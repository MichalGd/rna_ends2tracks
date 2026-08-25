from __future__ import annotations

import shutil
from collections import defaultdict
from pathlib import Path
from typing import Any

from .config import RunPlan, signature_for
from .external import event, require_tools, run
from .paths import workflow_asset
from .receipts import receipt_valid, write_receipt


def _threads(project: dict[str, Any]) -> int:
    return max(1, int(project.get("resources", {}).get("threads", 8)))


def preprocess(plan: RunPlan, results: Path, dry_run: bool = False, force: bool = False) -> None:
    module_dir = results / "02_alignment"
    qc_dir = results / "01_qc"
    log_dir = results / "provenance" / "logs"
    module_dir.mkdir(parents=True, exist_ok=True)
    qc_dir.mkdir(parents=True, exist_ok=True)
    inputs = [row["fastq_r1"] for row in plan.sample_rows]
    signature = signature_for(inputs, {"module": "preprocess", "project": plan.project, "reference": plan.reference,
                                       "lanes": plan.sample_rows})
    expected = [module_dir / sample["sample_id"] / f"{sample['sample_id']}.bam" for sample in plan.samples]
    expected += [path.with_suffix(".bam.bai") for path in expected]
    orientation_path = module_dir / "protocol_orientation.tsv"
    expected.append(orientation_path)
    if not force and receipt_valid(module_dir, signature):
        event(log_dir, "preprocess", "skipped", "Valid matching receipt")
        return
    if not dry_run:
        require_tools(["fastqc", "multiqc", "bbduk.sh", "STAR", "samtools"])
    threads = _threads(plan.project)
    trim_cfg = plan.project.get("preprocessing", {})
    trimq = str(trim_cfg.get("trim_quality", 10))
    minimum_length = str(trim_cfg.get("minimum_length", 20))
    adapters = str(trim_cfg.get("bbduk_reference", "adapters,polyA_T"))
    if adapters == "adapters,polyA_T":
        packaged = workflow_asset("resources/quantseq_rev_adapters.fa")
        adapters = str(packaged)

    lane_bams: dict[str, list[Path]] = defaultdict(list)
    orientation_rows: list[tuple[str, str, int, int, float, str]] = []
    raw_qc = qc_dir / "raw_fastqc"
    trimmed_qc = qc_dir / "trimmed_fastqc"
    trim_root = results / "01_qc" / "trimmed_fastq"
    raw_qc.mkdir(parents=True, exist_ok=True)
    trimmed_qc.mkdir(parents=True, exist_ok=True)
    trim_root.mkdir(parents=True, exist_ok=True)

    for lane in plan.sample_rows:
        sample = lane["sample_id"]
        lane_id = lane["lane_id"]
        token = f"{sample}.{lane_id}"
        lane_log = log_dir / "preprocess" / f"{token}.log"
        trimmed = trim_root / f"{token}.trimmed.fastq.gz"
        raw_lane_qc = raw_qc / token
        trimmed_lane_qc = trimmed_qc / token
        raw_lane_qc.mkdir(parents=True, exist_ok=True)
        trimmed_lane_qc.mkdir(parents=True, exist_ok=True)
        lane_prefix = module_dir / sample / "lanes" / token
        lane_prefix.parent.mkdir(parents=True, exist_ok=True)
        lane_bam = lane_prefix.with_suffix(".bam")
        lane_bams[sample].append(lane_bam)
        run(["fastqc", "--threads", str(threads), "--outdir", str(raw_lane_qc), lane["fastq_r1"]], lane_log, dry_run)
        run([
            "bbduk.sh", f"in={lane['fastq_r1']}", f"out={trimmed}", f"ref={adapters}",
            "int=f", "k=13", "ktrim=r", "useshortkmers=t", "mink=5", "qtrim=r",
            f"trimq={trimq}", f"minlength={minimum_length}", f"threads={threads}",
        ], lane_log, dry_run)
        run(["fastqc", "--threads", str(threads), "--outdir", str(trimmed_lane_qc), str(trimmed)], lane_log, dry_run)
        star_prefix = str(lane_prefix) + ".star."
        run([
            "STAR", "--runThreadN", str(threads), "--genomeDir", plan.reference["star_index"],
            "--readFilesIn", str(trimmed), "--readFilesCommand", "zcat",
            "--outFileNamePrefix", star_prefix, "--outSAMtype", "BAM", "Unsorted",
            "--outSAMattributes", "NH", "HI", "AS", "nM", "NM", "MD",
            "--outSAMattrRGline", f"ID:{lane_id}", f"SM:{sample}", f"LB:{sample}", "PL:ILLUMINA",
            "--outSAMstrandField", "intronMotif", "--quantMode", "GeneCounts",
        ], lane_log, dry_run)
        if not dry_run:
            count_path = Path(star_prefix + "ReadsPerGene.out.tab")
            forward = reverse = 0
            with count_path.open(encoding="utf-8") as handle:
                for line in handle:
                    fields = line.rstrip("\n").split("\t")
                    if len(fields) >= 4 and not fields[0].startswith("N_"):
                        forward += int(fields[2]); reverse += int(fields[3])
            informative = forward + reverse
            fraction = reverse / informative if informative else 0.0
            threshold = float(plan.project.get("protocol", {}).get("orientation_min_fraction", 0.75))
            status = "pass" if informative < 1000 or fraction >= threshold else "fail"
            orientation_rows.append((sample, lane_id, forward, reverse, fraction, status))
            if status == "fail":
                raise RuntimeError(
                    f"Protocol orientation mismatch for {token}: reverse-compatible fraction {fraction:.3f} < {threshold:.3f}"
                )
        run([
            "samtools", "sort", "-@", str(threads), "-o", str(lane_bam), star_prefix + "Aligned.out.bam"
        ], lane_log, dry_run)
        run(["samtools", "quickcheck", "-v", str(lane_bam)], lane_log, dry_run)

    for sample in plan.samples:
        sample_id = sample["sample_id"]
        sample_dir = module_dir / sample_id
        sample_bam = sample_dir / f"{sample_id}.bam"
        sample_log = log_dir / "preprocess" / f"{sample_id}.merge.log"
        bams = lane_bams[sample_id]
        if len(bams) == 1:
            if dry_run:
                run(["cp", str(bams[0]), str(sample_bam)], sample_log, True)
            else:
                temporary_bam = sample_dir / f".{sample_id}.bam.tmp"
                shutil.copy2(bams[0], temporary_bam)
                temporary_bam.replace(sample_bam)
        else:
            run(["samtools", "merge", "-@", str(threads), "-f", str(sample_bam), *map(str, bams)], sample_log, dry_run)
        run(["samtools", "index", "-@", str(threads), str(sample_bam)], sample_log, dry_run)
        if dry_run:
            run(["samtools", "flagstat", "-@", str(threads), "-O", "tsv", str(sample_bam)], sample_log, True)
        else:
            with (sample_dir / "flagstat.tsv").open("w", encoding="utf-8") as output:
                import subprocess
                subprocess.run(
                    ["samtools", "flagstat", "-@", str(threads), "-O", "tsv", str(sample_bam)],
                    stdout=output, stderr=subprocess.STDOUT, check=True, text=True,
                )
    run(["multiqc", "--force", "--outdir", str(qc_dir / "multiqc"), str(qc_dir), str(module_dir)],
        log_dir / "preprocess" / "multiqc.log", dry_run)
    if not dry_run:
        with orientation_path.open("w", encoding="utf-8", newline="") as handle:
            import csv
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["sample_id", "lane_id", "star_forward_count", "star_reverse_count", "reverse_compatible_fraction", "status"])
            writer.writerows(orientation_rows)
        write_receipt("preprocess", module_dir, signature, expected, ["rna-ends2tracks", "preprocess"])
    event(log_dir, "preprocess", "dry_run" if dry_run else "completed", f"Processed {len(plan.samples)} samples")
