from __future__ import annotations

import csv
import subprocess
from pathlib import Path

from .config import RunPlan, signature_for
from .external import event, require_tools, run
from .receipts import receipt_valid, write_receipt


def gene_expression(plan: RunPlan, results: Path, script_root: Path, dry_run: bool = False, force: bool = False) -> None:
    module_dir = results / "03_gene_expression"
    count_dir = module_dir / "counts"
    log_dir = results / "provenance" / "logs"
    module_dir.mkdir(parents=True, exist_ok=True)
    count_dir.mkdir(parents=True, exist_ok=True)
    bams = [results / "02_alignment" / s["sample_id"] / f"{s['sample_id']}.bam" for s in plan.samples]
    if not dry_run:
        missing = [str(path) for path in bams if not path.is_file()]
        if missing:
            raise RuntimeError("DGE requires completed sample BAMs: " + ", ".join(missing))
    signature = signature_for([*bams, plan.reference["gtf"]], {
        "module": "dge", "design": plan.project["design"], "samples": plan.samples, "contrasts": plan.contrasts,
        "reporting": plan.project.get("reporting", {}),
    }) if not dry_run else "dry-run"
    counts = count_dir / "gene_counts.featureCounts.tsv"
    result_index = module_dir / "deseq2" / "result_index.tsv"
    if not force and not dry_run and receipt_valid(module_dir, signature):
        event(log_dir, "dge", "skipped", "Valid matching receipt")
        return
    if not dry_run:
        require_tools(["featureCounts", "Rscript"])
    threads = max(1, int(plan.project.get("resources", {}).get("threads", 8)))
    fc_log = log_dir / "dge" / "featureCounts.log"
    run([
        "featureCounts", "-T", str(threads), "-a", plan.reference["gtf"], "-o", str(counts),
        "-t", "exon", "-g", "gene_id", "-s", "2", "--primary", *map(str, bams),
    ], fc_log, dry_run)
    r_script = script_root / "R" / "deseq2_all_pairs.R"
    run([
        "Rscript", str(r_script), "--counts", str(counts),
        "--samples", str(results / "00_metadata" / "validated_samples.tsv"),
        "--contrasts", str(results / "00_metadata" / "contrasts.tsv"),
        "--design", str(plan.project["design"]), "--outdir", str(module_dir / "deseq2"),
        "--fdr", str(plan.project.get("reporting", {}).get("fdr", 0.05)),
    ], log_dir / "dge" / "deseq2.log", dry_run)
    if not dry_run:
        write_receipt("dge", module_dir, signature, [counts, result_index], ["rna-ends2tracks", "dge"])
    event(log_dir, "dge", "dry_run" if dry_run else "completed", "Independent gene-expression branch")
