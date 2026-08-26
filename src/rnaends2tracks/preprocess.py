from __future__ import annotations

import csv
import json
import os
import shutil
import stat
import subprocess
import tempfile
from collections import defaultdict
from pathlib import Path
from typing import Any

from .config import RunPlan, signature_for
from .execution import run_bounded
from .external import event, require_tools, run
from .paths import workflow_asset
from .receipts import receipt_valid, write_receipt


def _fastqc_report(directory: Path, fastq: str | Path) -> Path:
    name = Path(fastq).name
    for suffix in (".gz", ".fastq", ".fq"):
        name = name.removesuffix(suffix)
    return directory / f"{name}_fastqc.html"


def _temporary_environment(project: dict[str, Any]) -> dict[str, str] | None:
    configured = str(project["resources"].get("temporary_directory", "") or "")
    if not configured:
        return None
    path = Path(configured).expanduser().resolve()
    path.mkdir(parents=True, exist_ok=True)
    return {"TMPDIR": str(path)}


def _remove_owned_temporary_tree(path: Path) -> None:
    """Remove a workflow-owned tree containing read-only copied package assets."""
    if not path.is_dir():
        return
    for root, directories, _files in os.walk(path):
        root_path = Path(root)
        root_path.chmod(root_path.stat().st_mode | stat.S_IWUSR | stat.S_IXUSR)
        for directory in directories:
            child = root_path / directory
            child.chmod(child.stat().st_mode | stat.S_IWUSR | stat.S_IXUSR)
    shutil.rmtree(path)


def preprocess(plan: RunPlan, results: Path, dry_run: bool = False, force: bool = False) -> None:
    """QC, trim, align and create the C0 unique-primary BAM universe."""
    module_dir = results / "02_alignment"
    qc_dir = results / "01_qc"
    log_dir = results / "logs"
    timing_dir = results / ".checkpoints" / "timings" / "preprocess"
    module_dir.mkdir(parents=True, exist_ok=True)
    qc_dir.mkdir(parents=True, exist_ok=True)
    reference_inputs: list[Path] = []
    if not dry_run:
        for reference in (plan.references or {plan.reference["assembly"]: plan.reference}).values():
            reference_inputs.extend(Path(reference[key]) for key in ("fasta", "gtf", "chrom_sizes"))
            star_index = Path(reference["star_index"])
            reference_inputs.extend(star_index / name for name in (
                "Genome", "SA", "SAindex", "chrName.txt", "chrLength.txt", "genomeParameters.txt",
            ))
    signature = signature_for([*[row["fastq_r1"] for row in plan.sample_rows], *reference_inputs], {
        "module": "preprocess", "protocol": plan.project.get("protocol", {}),
        "preprocessing": plan.project.get("preprocessing", {}),
        "references": plan.references or {plan.reference["assembly"]: plan.reference},
        "lanes": plan.sample_rows,
    }) if not dry_run else "dry-run"
    expected = [module_dir / sample["sample_id"] / f"{sample['sample_id']}.bam" for sample in plan.samples]
    expected += [path.with_suffix(".bam.bai") for path in expected]
    orientation_path = module_dir / "protocol_orientation.tsv"
    expected.extend([orientation_path, qc_dir / "multiqc" / "multiqc_report.html"])
    if not force and receipt_valid(module_dir, signature):
        event(log_dir, "preprocess", "skipped", "Valid matching module receipt")
        return
    if not dry_run:
        require_tools(["fastqc", "multiqc", "bbduk.sh", "STAR", "samtools"])

    resource = plan.project["resources"]["preprocess"]
    trim_cfg = plan.project.get("preprocessing", {})
    trimq = str(trim_cfg.get("trim_quality", 10))
    minimum_length = str(trim_cfg.get("minimum_length", 20))
    adapters = str(trim_cfg.get("bbduk_reference", "adapters,polyA_T"))
    if adapters == "adapters,polyA_T":
        adapters = str(workflow_asset("resources/quantseq_rev_adapters.fa"))
    tool_env = _temporary_environment(plan.project)
    raw_qc = qc_dir / "raw_fastqc"
    trimmed_qc = qc_dir / "trimmed_fastqc"
    trim_root = qc_dir / "trimmed_fastq"
    for directory in (raw_qc, trimmed_qc, trim_root):
        directory.mkdir(parents=True, exist_ok=True)

    lane_bams: dict[str, list[Path]] = defaultdict(list)
    contexts: list[dict[str, Any]] = []
    trim_jobs: list[tuple[str, Any]] = []
    for lane in plan.sample_rows:
        sample = lane["sample_id"]
        token = f"{sample}.{lane['technical_replicate_id']}.{lane['lane_id']}"
        lane_prefix = module_dir / sample / "lanes" / token
        lane_prefix.parent.mkdir(parents=True, exist_ok=True)
        lane_bam = lane_prefix.parent / f"{token}.bam"
        trimmed = trim_root / f"{token}.trimmed.fastq.gz"
        lane_bams[sample].append(lane_bam)
        context = {"lane": lane, "token": token, "lane_prefix": lane_prefix,
                   "lane_bam": lane_bam, "trimmed": trimmed}
        contexts.append(context)

        def trim_worker(context=context):
            lane = context["lane"]
            token = context["token"]
            trimmed = context["trimmed"]
            lane_log = log_dir / "preprocess" / f"{token}.trim.log"
            receipt_dir = trim_root / ".receipts" / token
            trim_signature = signature_for([lane["fastq_r1"]], {
                "module": "qc_trim_lane", "lane": lane, "preprocessing": trim_cfg,
                "resources": {key: resource[key] for key in (
                    "fastqc_threads", "bbduk_threads", "bbduk_memory_gb")},
            }) if not dry_run else "dry-run"
            if not force and not dry_run and receipt_valid(receipt_dir, trim_signature):
                return trimmed
            raw_lane_qc = raw_qc / token
            trimmed_lane_qc = trimmed_qc / token
            raw_lane_qc.mkdir(parents=True, exist_ok=True)
            trimmed_lane_qc.mkdir(parents=True, exist_ok=True)
            run(["fastqc", "--threads", str(resource["fastqc_threads"]), "--outdir", str(raw_lane_qc),
                 lane["fastq_r1"]], lane_log, dry_run, env=tool_env)
            temporary = trimmed.with_name(f".{trimmed.name}.tmp.fastq.gz")
            run(["bbduk.sh", f"-Xmx{resource['bbduk_memory_gb']}g", f"in={lane['fastq_r1']}",
                 f"out={temporary}", f"ref={adapters}", "int=f", "k=13", "ktrim=r",
                 "useshortkmers=t", "mink=5", "qtrim=r", f"trimq={trimq}",
                 f"minlength={minimum_length}", f"threads={resource['bbduk_threads']}"],
                lane_log, dry_run, env=tool_env)
            if not dry_run:
                temporary.replace(trimmed)
            run(["fastqc", "--threads", str(resource["fastqc_threads"]), "--outdir", str(trimmed_lane_qc),
                 str(trimmed)], lane_log, dry_run, env=tool_env)
            if not dry_run:
                write_receipt(
                    "qc_trim_lane", receipt_dir, trim_signature,
                    [trimmed, _fastqc_report(raw_lane_qc, lane["fastq_r1"]),
                     _fastqc_report(trimmed_lane_qc, trimmed)],
                    ["rna-ends2tracks", "qc_trim", token],
                )
            return trimmed

        trim_jobs.append((token, trim_worker))

    run_bounded("qc_and_trim", trim_jobs, resource["trim_parallel_jobs"], timing_dir / "trim")

    alignment_jobs: list[tuple[str, Any]] = []
    for context in contexts:
        def alignment_worker(context=context):
            lane = context["lane"]
            token = context["token"]
            lane_prefix = context["lane_prefix"]
            lane_bam = context["lane_bam"]
            trimmed = context["trimmed"]
            reference = plan.reference_for(lane["genome"])
            lane_log = log_dir / "preprocess" / f"{token}.star.log"
            receipt_dir = module_dir / lane["sample_id"] / "lanes" / ".receipts" / token
            orientation_json = receipt_dir / "orientation.json"
            star_index = Path(reference["star_index"])
            lane_signature = signature_for([lane["fastq_r1"], *[
                star_index / name for name in (
                    "Genome", "SA", "SAindex", "chrName.txt", "chrLength.txt", "genomeParameters.txt",
                )
            ]], {
                "module": "star_lane", "lane": lane, "protocol": plan.project.get("protocol", {}),
                "preprocessing": trim_cfg, "reference": reference,
                "resources": {key: resource[key] for key in (
                    "star_threads", "star_memory_gb", "samtools_threads",
                    "samtools_sort_memory_per_thread_gb")},
            }) if not dry_run else "dry-run"
            if not force and not dry_run and receipt_valid(receipt_dir, lane_signature):
                return json.loads(orientation_json.read_text(encoding="utf-8"))
            star_prefix = str(lane_prefix) + ".star."
            run(["STAR", "--runThreadN", str(resource["star_threads"]),
                 "--limitGenomeGenerateRAM", str(resource["star_memory_gb"] * 1024**3),
                 "--genomeDir", reference["star_index"], "--readFilesIn", str(trimmed),
                 "--readFilesCommand", "zcat", "--outFileNamePrefix", star_prefix,
                 "--outSAMtype", "BAM", "Unsorted", "--outSAMattributes", "NH", "HI", "AS", "nM", "NM", "MD",
                 "--outSAMattrRGline", f"ID:{token}", f"SM:{lane['sample_id']}",
                 f"LB:{lane['sample_id']}.{lane['technical_replicate_id']}", "PL:ILLUMINA",
                 "--outSAMstrandField", "intronMotif", "--quantMode", "GeneCounts"],
                lane_log, dry_run, env=tool_env)
            orientation = {"sample_id": lane["sample_id"],
                           "technical_replicate_id": lane["technical_replicate_id"],
                           "lane_id": lane["lane_id"], "star_forward_count": 0,
                           "star_reverse_count": 0, "reverse_compatible_fraction": 0.0,
                           "status": "dry_run" if dry_run else "pass"}
            if not dry_run:
                forward = reverse = 0
                with Path(star_prefix + "ReadsPerGene.out.tab").open(encoding="utf-8") as handle:
                    for line in handle:
                        fields = line.rstrip("\n").split("\t")
                        if len(fields) >= 4 and not fields[0].startswith("N_"):
                            forward += int(fields[2])
                            reverse += int(fields[3])
                informative = forward + reverse
                fraction = reverse / informative if informative else 0.0
                threshold = float(plan.project.get("protocol", {}).get("orientation_min_fraction", 0.75))
                status = "pass" if informative < 1000 or fraction >= threshold else "fail"
                orientation.update({"star_forward_count": forward, "star_reverse_count": reverse,
                                    "reverse_compatible_fraction": fraction, "status": status})
                if status == "fail":
                    raise RuntimeError(f"Protocol orientation mismatch for {token}: reverse-compatible fraction "
                                       f"{fraction:.3f} < {threshold:.3f}")
            temporary_bam = lane_bam.with_name(f".{lane_bam.name}.tmp")
            run(["samtools", "sort", "-@", str(resource["samtools_threads"]),
                 "-m", f"{resource['samtools_sort_memory_per_thread_gb']}G", "-o", str(temporary_bam),
                 star_prefix + "Aligned.out.bam"], lane_log, dry_run, env=tool_env)
            if not dry_run:
                temporary_bam.replace(lane_bam)
            run(["samtools", "quickcheck", "-v", str(lane_bam)], lane_log, dry_run, env=tool_env)
            if not dry_run:
                receipt_dir.mkdir(parents=True, exist_ok=True)
                temporary_orientation = receipt_dir / ".orientation.json.tmp"
                temporary_orientation.write_text(json.dumps(orientation, indent=2, sort_keys=True) + "\n",
                                                 encoding="utf-8")
                temporary_orientation.replace(orientation_json)
                write_receipt("star_lane", receipt_dir, lane_signature, [lane_bam, orientation_json],
                              ["rna-ends2tracks", "alignment", token])
            return orientation

        alignment_jobs.append((context["token"], alignment_worker))

    orientation_rows = run_bounded("star_and_sort", alignment_jobs, resource["star_parallel_jobs"],
                                   timing_dir / "alignment")

    merge_jobs: list[tuple[str, Any]] = []
    for sample in plan.samples:
        sample_id = sample["sample_id"]
        sample_dir = module_dir / sample_id
        sample_dir.mkdir(parents=True, exist_ok=True)
        sample_bam = sample_dir / f"{sample_id}.bam"
        bams = lane_bams[sample_id]

        def merge_worker(sample_id=sample_id, sample_dir=sample_dir, sample_bam=sample_bam, bams=bams):
            sample_log = log_dir / "preprocess" / f"{sample_id}.merge.log"
            receipt_dir = sample_dir / ".merge_receipt"
            sample_signature = signature_for(bams, {"module": "preprocess_merge", "sample_id": sample_id}) if not dry_run else "dry-run"
            flagstat = sample_dir / "flagstat.tsv"
            if not force and not dry_run and receipt_valid(receipt_dir, sample_signature):
                return sample_bam
            temporary_bam = sample_dir / f".{sample_id}.bam.tmp"
            merged_all = sample_dir / f"{sample_id}.all_alignments.bam"
            merge_source = bams[0]
            if len(bams) > 1:
                temporary_merged = sample_dir / f".{sample_id}.all_alignments.bam.tmp"
                run(["samtools", "merge", "-@", str(resource["samtools_threads"]), "-f",
                     str(temporary_merged), *map(str, bams)], sample_log, dry_run, env=tool_env)
                if not dry_run:
                    temporary_merged.replace(merged_all)
                merge_source = merged_all
            # C0: mapped primary NH=1. Duplicate flags are deliberately retained.
            run(["samtools", "view", "-@", str(resource["samtools_threads"]), "-b", "-F", "2308",
                 "-e", "[NH] == 1", "-o", str(temporary_bam), str(merge_source)],
                sample_log, dry_run, env=tool_env)
            if not dry_run:
                temporary_bam.replace(sample_bam)
            temporary_index = sample_dir / f".{sample_id}.bam.bai.tmp"
            run(["samtools", "index", "-@", str(resource["samtools_threads"]), "-o", str(temporary_index),
                 str(sample_bam)], sample_log, dry_run, env=tool_env)
            if not dry_run:
                temporary_index.replace(sample_bam.with_suffix(".bam.bai"))
                temporary_flagstat = sample_dir / ".flagstat.tsv.tmp"
                with temporary_flagstat.open("w", encoding="utf-8") as output:
                    subprocess.run(["samtools", "flagstat", "-@", str(resource["samtools_threads"]), "-O", "tsv",
                                    str(sample_bam)], stdout=output, stderr=subprocess.STDOUT, check=True, text=True,
                                   env=None if tool_env is None else {**os.environ, **tool_env})
                temporary_flagstat.replace(flagstat)
                write_receipt("preprocess_merge", receipt_dir, sample_signature,
                              [sample_bam, sample_bam.with_suffix(".bam.bai"), flagstat],
                              ["rna-ends2tracks", "preprocess_merge", sample_id])
            else:
                run(["samtools", "flagstat", "-@", str(resource["samtools_threads"]), "-O", "tsv", str(sample_bam)],
                    sample_log, True, env=tool_env)
            return sample_bam

        merge_jobs.append((sample_id, merge_worker))

    run_bounded("preprocess_merges", merge_jobs, resource["merge_parallel_jobs"], timing_dir / "merges")
    # Frozen shared environments make MultiQC's copied template directories read-only.
    # Disable MultiQC's own cleanup and remove only our dedicated temporary tree.
    multiqc_temp_parent = results / ".checkpoints" / "multiqc_tmp"
    multiqc_temp_parent.mkdir(parents=True, exist_ok=True)
    if dry_run:
        multiqc_temp = multiqc_temp_parent / "dry-run"
    else:
        multiqc_temp = Path(tempfile.mkdtemp(prefix="run.", dir=multiqc_temp_parent))
    multiqc_env = {**(tool_env or {}), "TMPDIR": str(multiqc_temp)}
    run(["multiqc", "--force", "--no-clean-up", "--outdir", str(qc_dir / "multiqc"),
         str(qc_dir), str(module_dir)],
        log_dir / "preprocess" / "multiqc.log", dry_run, env=multiqc_env)
    if not dry_run:
        _remove_owned_temporary_tree(multiqc_temp)
    if not dry_run:
        temporary_orientation = orientation_path.with_name(f".{orientation_path.name}.tmp")
        with temporary_orientation.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=["sample_id", "technical_replicate_id", "lane_id",
                "star_forward_count", "star_reverse_count", "reverse_compatible_fraction", "status"],
                delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(orientation_rows)
        temporary_orientation.replace(orientation_path)
        write_receipt("preprocess", module_dir, signature, expected, ["rna-ends2tracks", "preprocess"])
    event(log_dir, "preprocess", "dry_run" if dry_run else "completed",
          f"Processed {len(plan.samples)} samples from {len(plan.sample_rows)} lanes; "
          f"trim_parallel_jobs={resource['trim_parallel_jobs']}; star_parallel_jobs={resource['star_parallel_jobs']}")
