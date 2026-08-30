from __future__ import annotations

import csv
import json
import os
import shutil
import stat
import tempfile
from collections import defaultdict
from concurrent.futures import Future, ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Any

from .config import RunPlan, signature_for
from .execution import run_bounded
from .external import event, progress_events, require_tools, run, run_to_path
from .paths import workflow_asset
from .receipts import receipt_valid, write_receipt
from .tracks import c0_tracks_enabled, make_c0_tracks_for_sample


def _fastqc_report(directory: Path, fastq: str | Path) -> Path:
    name = Path(fastq).name
    for suffix in (".gz", ".fastq", ".fq"):
        name = name.removesuffix(suffix)
    return directory / f"{name}_fastqc.html"


def _fastq_inputs(lane: dict[str, str]) -> list[str]:
    return [lane["fastq_r1"], *([lane["fastq_r2"]] if lane.get("fastq_r2") else [])]


def _fastq_screen_mate(report: Path, lane: dict[str, str]) -> str:
    if report.parent.name in {"R1", "R2"}:
        return report.parent.name
    for mate, key in (("R1", "fastq_r1"), ("R2", "fastq_r2")):
        value = lane.get(key, "")
        if not value:
            continue
        name = Path(value).name
        for suffix in (".fastq.gz", ".fq.gz", ".fastq", ".fq", ".gz"):
            if name.endswith(suffix):
                name = name[:-len(suffix)]
                break
        if report.name.startswith(f"{name}_screen"):
            return mate
    return "UNKNOWN"


def _fastq_screen_metrics(
    report: Path, lane: dict[str, str], token: str,
) -> list[dict[str, str]]:
    lines = report.read_text(encoding="utf-8", errors="replace").splitlines()
    header = next((index for index, line in enumerate(lines) if line.startswith("Genome\t")), None)
    if header is None:
        raise RuntimeError(f"FastQ Screen report has no tabular Genome header: {report}")
    reader = csv.DictReader(lines[header:], delimiter="\t")
    rows: list[dict[str, str]] = []
    for value in reader:
        if not value.get("Genome"):
            continue
        rows.append({
            "sample_id": lane["sample_id"],
            "technical_replicate_id": lane["technical_replicate_id"],
            "lane_id": lane["lane_id"], "lane_token": token,
            "mate": _fastq_screen_mate(report, lane), "database": value["Genome"],
            "reads_processed": value.get("#Reads_processed", ""),
            "pct_unmapped": value.get("%Unmapped", ""),
            "pct_one_hit_one_library": value.get("%One_hit_one_library", ""),
            "pct_multiple_hits_one_library": value.get("%Multiple_hits_one_library", ""),
            "pct_one_hit_multiple_libraries": value.get("%One_hit_multiple_libraries", ""),
            "pct_multiple_hits_multiple_libraries": value.get("%Multiple_hits_multiple_libraries", ""),
            "source_report": str(report),
        })
    if not rows:
        raise RuntimeError(f"FastQ Screen report contains no database rows: {report}")
    return rows


def _run_fastq_screen(
    plan: RunPlan, results: Path, dry_run: bool, force: bool,
    tool_env: dict[str, str] | None,
) -> Path:
    """Run lane/mate-aware FastQ Screen and publish a compact status table."""
    settings = plan.project.get("preprocessing", {}).get("fastq_screen", {})
    root = results / "01_qc" / "fastq_screen"
    summary = root / "fastq_screen_summary.tsv"
    metrics_path = root / "fastq_screen_metrics.tsv"
    root.mkdir(parents=True, exist_ok=True)
    config = Path(str(settings.get("config", ""))) if settings.get("config") else None
    enabled = bool(settings.get("enabled", True))
    missing_action = str(settings.get("missing_action", "warn"))
    rows: list[dict[str, str]] = []
    if not enabled or config is None or not config.is_file():
        if enabled and missing_action == "error":
            raise RuntimeError(
                "FastQ Screen is enabled but FASTQ_SCREEN_CONFIG does not identify a readable database config"
            )
        status = "DISABLED" if not enabled else "SKIPPED_MISSING_CONFIG"
        for lane in plan.sample_rows:
            rows.append({
                "sample_id": lane["sample_id"], "technical_replicate_id": lane["technical_replicate_id"],
                "lane_id": lane["lane_id"], "layout": lane.get("library_layout", "SE"),
                "mates": "R1,R2" if lane.get("fastq_r2") else "R1", "status": status,
                "output_directory": "", "config": str(config or ""),
            })
        if enabled:
            event(results / "logs", "fastq_screen", "warning",
                  "FASTQ_SCREEN_CONFIG is unavailable; contamination screening was explicitly skipped")
    else:
        if not dry_run:
            require_tools(["fastq_screen"])
        resource = plan.project["resources"]["preprocess"]
        jobs: list[tuple[str, Any]] = []
        for lane in plan.sample_rows:
            token = f"{lane['sample_id']}.{lane['technical_replicate_id']}.{lane['lane_id']}"
            outdir = root / token

            def worker(lane=lane, token=token, outdir=outdir):
                inputs = _fastq_inputs(lane)
                receipt_dir = outdir / ".receipt"
                marker = outdir / "fastq_screen.complete.json"
                signature = signature_for([*inputs, config], {
                    "module": "fastq_screen", "subset": settings.get("subset", 200000),
                    "threads": resource["fastq_screen_threads"], "layout": lane.get("library_layout", "SE"),
                }) if not dry_run else "dry-run"
                if not force and not dry_run and receipt_valid(receipt_dir, signature):
                    reports = sorted(outdir.glob("*/*_screen.txt"))
                    return {
                        "sample_id": lane["sample_id"], "technical_replicate_id": lane["technical_replicate_id"],
                        "lane_id": lane["lane_id"], "layout": lane.get("library_layout", "SE"),
                        "mates": "R1,R2" if lane.get("fastq_r2") else "R1", "status": "PASS",
                        "output_directory": str(outdir), "config": str(config),
                        "text_reports": ",".join(map(str, reports)),
                        "_metrics": [item for report in reports
                                     for item in _fastq_screen_metrics(report, lane, token)],
                    }
                outdir.mkdir(parents=True, exist_ok=True)
                commands: list[list[str]] = []
                for mate, fastq in zip(("R1", "R2"), inputs):
                    mate_dir = outdir / mate
                    mate_dir.mkdir(parents=True, exist_ok=True)
                    command = [
                        "fastq_screen", "--conf", str(config), "--outdir", str(mate_dir),
                        "--threads", str(resource["fastq_screen_threads"]),
                        "--subset", str(settings.get("subset", 200000)), "--force", fastq,
                    ]
                    commands.append(command)
                    run(command, results / "logs" / "preprocess" / f"{token}.{mate}.fastq_screen.log",
                        dry_run, env=tool_env)
                if not dry_run:
                    reports = sorted(outdir.glob("*/*_screen.txt"))
                    if len(reports) != len(inputs):
                        raise RuntimeError(
                            f"FastQ Screen produced {len(reports)} text reports for {len(inputs)} inputs: {token}"
                        )
                    marker.write_text(json.dumps({"status": "PASS", "inputs": inputs,
                                                  "commands": commands}, indent=2) + "\n",
                                      encoding="utf-8")
                    write_receipt(
                        "fastq_screen_lane", receipt_dir, signature, [marker, *reports],
                        ["rna-ends2tracks", "fastq_screen", token, *inputs],
                    )
                else:
                    reports = []
                return {
                    "sample_id": lane["sample_id"], "technical_replicate_id": lane["technical_replicate_id"],
                    "lane_id": lane["lane_id"], "layout": lane.get("library_layout", "SE"),
                    "mates": "R1,R2" if lane.get("fastq_r2") else "R1",
                    "status": "DRY_RUN" if dry_run else "PASS", "output_directory": str(outdir),
                    "config": str(config), "text_reports": ",".join(map(str, reports)),
                    "_metrics": [item for report in reports
                                 for item in _fastq_screen_metrics(report, lane, token)],
                }

            jobs.append((token, worker))
        rows = run_bounded(
            "fastq_screen", jobs, resource["fastq_screen_parallel_jobs"],
            results / ".checkpoints" / "timings" / "preprocess" / "fastq_screen",
            progress=progress_events(results / "logs", "fastq_screen", len(jobs), "lane"),
        )
    metrics: list[dict[str, str]] = []
    for row in rows:
        metrics.extend(row.pop("_metrics", []))
    temporary = summary.with_name(f".{summary.name}.tmp")
    with temporary.open("w", encoding="utf-8", newline="") as handle:
        fields = ["sample_id", "technical_replicate_id", "lane_id", "layout", "mates", "status",
                  "output_directory", "text_reports", "config"]
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)
    temporary.replace(summary)
    metric_fields = [
        "sample_id", "technical_replicate_id", "lane_id", "lane_token", "mate", "database",
        "reads_processed", "pct_unmapped", "pct_one_hit_one_library",
        "pct_multiple_hits_one_library", "pct_one_hit_multiple_libraries",
        "pct_multiple_hits_multiple_libraries", "source_report",
    ]
    metrics_temporary = metrics_path.with_name(f".{metrics_path.name}.tmp")
    with metrics_temporary.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=metric_fields, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(metrics)
    metrics_temporary.replace(metrics_path)
    return summary


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


def _c0_overlap_workers(plan: RunPlan) -> int:
    """Return a conservative track-worker count that fits beside sample merging."""
    if not c0_tracks_enabled(plan):
        return 0
    resources = plan.project["resources"]
    preprocess_resources = resources["preprocess"]
    track_resources = resources["tracks"]
    merge_workers = min(preprocess_resources["merge_parallel_jobs"], len(plan.samples))
    remaining_threads = (
        resources["total_threads"]
        - merge_workers * preprocess_resources["samtools_threads"]
    )
    remaining_memory = (
        resources["total_memory_gb"]
        - merge_workers * preprocess_resources["merge_memory_gb"]
    )
    by_threads = remaining_threads // track_resources["samtools_threads"]
    by_memory = remaining_memory // track_resources["memory_gb"]
    return max(0, min(
        track_resources["parallel_jobs"], len(plan.samples), by_threads, by_memory,
    ))


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
    fastq_inputs = [path for row in plan.sample_rows for path in _fastq_inputs(row)]
    signature = signature_for([*fastq_inputs, *reference_inputs], {
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
    fastq_screen_summary = _run_fastq_screen(plan, results, dry_run, force, tool_env)
    if not dry_run:
        expected.extend([fastq_screen_summary, fastq_screen_summary.with_name("fastq_screen_metrics.tsv")])

    lane_bams: dict[str, list[Path]] = defaultdict(list)
    contexts: list[dict[str, Any]] = []
    trim_jobs: list[tuple[str, Any]] = []
    for lane in plan.sample_rows:
        sample = lane["sample_id"]
        token = f"{sample}.{lane['technical_replicate_id']}.{lane['lane_id']}"
        lane_prefix = module_dir / sample / "lanes" / token
        lane_prefix.parent.mkdir(parents=True, exist_ok=True)
        lane_bam = lane_prefix.parent / f"{token}.bam"
        paired = lane.get("library_layout", "SE") == "PE"
        trimmed_r1 = trim_root / (f"{token}.R1.trimmed.fastq.gz" if paired else f"{token}.trimmed.fastq.gz")
        trimmed_r2 = trim_root / f"{token}.R2.trimmed.fastq.gz" if paired else None
        lane_bams[sample].append(lane_bam)
        context = {"lane": lane, "token": token, "lane_prefix": lane_prefix,
                   "lane_bam": lane_bam, "trimmed_r1": trimmed_r1, "trimmed_r2": trimmed_r2}
        contexts.append(context)

        def trim_worker(context=context):
            lane = context["lane"]
            token = context["token"]
            trimmed_r1 = context["trimmed_r1"]
            trimmed_r2 = context["trimmed_r2"]
            lane_log = log_dir / "preprocess" / f"{token}.trim.log"
            receipt_dir = trim_root / ".receipts" / token
            trim_signature = signature_for(_fastq_inputs(lane), {
                "module": "qc_trim_lane", "lane": lane, "preprocessing": trim_cfg,
                "resources": {key: resource[key] for key in (
                    "fastqc_threads", "bbduk_threads", "bbduk_memory_gb")},
            }) if not dry_run else "dry-run"
            if not force and not dry_run and receipt_valid(receipt_dir, trim_signature):
                return trimmed_r1
            raw_lane_qc = raw_qc / token
            trimmed_lane_qc = trimmed_qc / token
            raw_lane_qc.mkdir(parents=True, exist_ok=True)
            trimmed_lane_qc.mkdir(parents=True, exist_ok=True)
            run(["fastqc", "--threads", str(resource["fastqc_threads"]), "--outdir", str(raw_lane_qc),
                 *_fastq_inputs(lane)], lane_log, dry_run, env=tool_env)
            temporary_r1 = trimmed_r1.with_name(f".{trimmed_r1.name}.tmp.fastq.gz")
            temporary_r2 = (trimmed_r2.with_name(f".{trimmed_r2.name}.tmp.fastq.gz")
                            if trimmed_r2 is not None else None)
            adapter_r1 = (trimmed_r1.with_name(f".{trimmed_r1.name}.adapter.tmp.fastq.gz")
                          if trimmed_r2 is not None else temporary_r1)
            adapter_r2 = (trimmed_r2.with_name(f".{trimmed_r2.name}.adapter.tmp.fastq.gz")
                          if trimmed_r2 is not None else None)
            io_args = ([f"in1={lane['fastq_r1']}", f"in2={lane['fastq_r2']}",
                        f"out1={adapter_r1}", f"out2={adapter_r2}"]
                       if trimmed_r2 is not None else [f"in={lane['fastq_r1']}", f"out={temporary_r1}"])
            run(["bbduk.sh", f"-Xmx{resource['bbduk_memory_gb']}g", *io_args,
                 f"ref={adapters}", "int=f", "k=13", "ktrim=r",
                 "useshortkmers=t", "mink=5", "qtrim=r", f"trimq={trimq}",
                 f"minlength={minimum_length}", f"threads={resource['bbduk_threads']}"],
                lane_log, dry_run, env=tool_env)
            r2_trim = int(trim_cfg.get("pe_r2_trim_5p", 12)) if trimmed_r2 is not None else 0
            if trimmed_r2 is not None and r2_trim:
                # BBTools documents skipr1 as applying to force trimming. This
                # synchronized pass therefore removes only the random-primer-
                # derived R2 prefix while preserving pairing and R1.
                run([
                    "bbduk.sh", f"-Xmx{resource['bbduk_memory_gb']}g",
                    f"in1={adapter_r1}", f"in2={adapter_r2}",
                    f"out1={temporary_r1}", f"out2={temporary_r2}",
                    f"ftl={r2_trim}", "skipr1=t", "int=f",
                    f"threads={resource['bbduk_threads']}",
                ], lane_log, dry_run, env=tool_env)
                if not dry_run:
                    adapter_r1.unlink(missing_ok=True)
                    adapter_r2.unlink(missing_ok=True)
            elif trimmed_r2 is not None and not dry_run:
                adapter_r1.replace(temporary_r1)
                adapter_r2.replace(temporary_r2)
            if not dry_run:
                temporary_r1.replace(trimmed_r1)
                if temporary_r2 is not None and trimmed_r2 is not None:
                    temporary_r2.replace(trimmed_r2)
            run(["fastqc", "--threads", str(resource["fastqc_threads"]), "--outdir", str(trimmed_lane_qc),
                 str(trimmed_r1), *([str(trimmed_r2)] if trimmed_r2 is not None else [])],
                lane_log, dry_run, env=tool_env)
            if not dry_run:
                outputs = [trimmed_r1, *([trimmed_r2] if trimmed_r2 is not None else []),
                           *[_fastqc_report(raw_lane_qc, path) for path in _fastq_inputs(lane)],
                           _fastqc_report(trimmed_lane_qc, trimmed_r1)]
                if trimmed_r2 is not None:
                    outputs.append(_fastqc_report(trimmed_lane_qc, trimmed_r2))
                write_receipt(
                    "qc_trim_lane", receipt_dir, trim_signature, outputs,
                    ["rna-ends2tracks", "qc_trim", token],
                )
            return trimmed_r1

        trim_jobs.append((token, trim_worker))

    run_bounded(
        "qc_and_trim", trim_jobs, resource["trim_parallel_jobs"], timing_dir / "trim",
        progress=progress_events(log_dir, "alignment", len(trim_jobs), "QC/trim lane"),
    )

    alignment_jobs: list[tuple[str, Any]] = []
    for context in contexts:
        def alignment_worker(context=context):
            lane = context["lane"]
            token = context["token"]
            lane_prefix = context["lane_prefix"]
            lane_bam = context["lane_bam"]
            trimmed_r1 = context["trimmed_r1"]
            trimmed_r2 = context["trimmed_r2"]
            reference = plan.reference_for(lane["genome"])
            lane_log = log_dir / "preprocess" / f"{token}.star.log"
            receipt_dir = module_dir / lane["sample_id"] / "lanes" / ".receipts" / token
            orientation_json = receipt_dir / "orientation.json"
            star_index = Path(reference["star_index"])
            lane_signature = signature_for([*_fastq_inputs(lane), *[
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
                 "--genomeDir", reference["star_index"], "--readFilesIn", str(trimmed_r1),
                 *([str(trimmed_r2)] if trimmed_r2 is not None else []),
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

    orientation_rows = run_bounded(
        "star_and_sort", alignment_jobs, resource["star_parallel_jobs"], timing_dir / "alignment",
        progress=progress_events(log_dir, "alignment", len(alignment_jobs), "STAR lane"),
    )

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
                run_to_path(
                    ["samtools", "flagstat", "-@", str(resource["samtools_threads"]), "-O", "tsv",
                     str(sample_bam)],
                    temporary_flagstat, sample_log, env=tool_env,
                )
                temporary_flagstat.replace(flagstat)
                write_receipt("preprocess_merge", receipt_dir, sample_signature,
                              [sample_bam, sample_bam.with_suffix(".bam.bai"), flagstat],
                              ["rna-ends2tracks", "preprocess_merge", sample_id])
            else:
                run(["samtools", "flagstat", "-@", str(resource["samtools_threads"]), "-O", "tsv", str(sample_bam)],
                    sample_log, True, env=tool_env)
            return sample_bam

        merge_jobs.append((sample_id, merge_worker))

    overlap_workers = _c0_overlap_workers(plan)
    early_c0_enabled = c0_tracks_enabled(plan)
    sample_by_id = {sample["sample_id"]: sample for sample in plan.samples}
    merge_progress = progress_events(log_dir, "alignment", len(merge_jobs), "C0 sample")
    if dry_run:
        run_bounded(
            "preprocess_merges", merge_jobs, resource["merge_parallel_jobs"], timing_dir / "merges",
            progress=merge_progress,
        )
        if early_c0_enabled:
            event(
                log_dir, "c0_tracks", "dry_run",
                f"Would stream sample-ready C0 tracks with {overlap_workers} overlap worker(s); "
                "zero means defer to the dedicated c0_tracks stage",
            )
    elif overlap_workers:
        event(
            log_dir, "c0_tracks", "started",
            f"Streaming sample-ready C0 tracks with {overlap_workers} worker(s) while BAM merging continues",
        )
        track_progress = progress_events(log_dir, "c0_tracks", len(merge_jobs), "sample")
        track_failures: list[tuple[str, Exception]] = []
        with ThreadPoolExecutor(
            max_workers=overlap_workers, thread_name_prefix="sample_ready_c0_tracks",
        ) as track_pool:
            track_futures: dict[Future[Any], str] = {}

            def submit_c0_tracks(sample_id: str, _sample_bam: Path) -> None:
                track_futures[
                    track_pool.submit(
                        make_c0_tracks_for_sample,
                        plan, results, sample_by_id[sample_id], False,
                    )
                ] = sample_id

            run_bounded(
                "preprocess_merges", merge_jobs, resource["merge_parallel_jobs"],
                timing_dir / "merges", progress=merge_progress,
                on_completed=submit_c0_tracks,
            )
            for future in as_completed(track_futures):
                sample_id = track_futures[future]
                try:
                    future.result()
                    track_progress(sample_id, "completed")
                except Exception as exc:  # noqa: BLE001 - report every sample failure
                    track_failures.append((sample_id, exc))
                    track_progress(sample_id, "failed")
        if track_failures:
            detail = "; ".join(f"{sample_id}: {exc}" for sample_id, exc in track_failures)
            raise RuntimeError(f"sample-ready C0 track worker failures: {detail}")
        event(
            log_dir, "c0_tracks", "completed",
            f"Published sample-ready C0 tracks for {len(track_futures)} sample(s)",
        )
    else:
        if early_c0_enabled:
            event(
                log_dir, "c0_tracks", "deferred",
                "Combined CPU/RAM budget leaves no safe overlap capacity; "
                "C0 tracks will run in the dedicated c0_tracks stage",
            )
        run_bounded(
            "preprocess_merges", merge_jobs, resource["merge_parallel_jobs"], timing_dir / "merges",
            progress=merge_progress,
        )
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
