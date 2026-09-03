from __future__ import annotations

import csv
from pathlib import Path

from .apa_mcell import _write_differential_pcpa
from .config import RunPlan, signature_for
from .external import event, require_tools
from .receipts import write_receipt
from .statistics import run_r_contrasts


def apa_a2(
    plan: RunPlan,
    results: Path,
    script_root: Path,
    dry_run: bool = False,
    force: bool = False,
) -> None:
    """Run corrected raw-count PAU effects in an independent APA-A2 stage."""
    log_dir = results / "logs"
    if dry_run:
        event(
            log_dir,
            "apa_a2",
            "dry_run",
            "Would run independent DEXSeq tests with raw-count within-gene PAU effects",
        )
        return
    require_tools(["Rscript"])
    active_root = results / "04_active_pas"
    output_root = results / "06b_apa_a2_corrected"
    outputs: list[Path] = []
    inputs: list[Path] = []
    reporting = plan.project["reporting"]

    for genome in plan.references or {plan.reference["assembly"]: plan.reference}:
        genome_samples = [sample for sample in plan.samples if sample["genome"] == genome]
        genome_rows = [row for row in plan.sample_rows if row["genome"] == genome]
        genome_contrasts = [contrast for contrast in plan.contrasts if contrast["genome"] == genome]
        if not genome_contrasts:
            continue
        counts = active_root / genome / "C3_active_pas_counts.tsv"
        catalog = active_root / genome / "active_pas_catalog.tsv"
        inputs.extend([counts, catalog])
        genome_plan = RunPlan(
            plan.project,
            genome_samples,
            genome_rows,
            genome_contrasts,
            plan.reference_for(genome),
            {genome: plan.reference_for(genome)},
        )
        outdir = output_root / genome / "dexseq_a2"
        index_path = outdir / "result_index.tsv"
        run_r_contrasts(
            module=f"apa_a2_{genome}",
            plan=genome_plan,
            results=results,
            script=script_root / "R" / "dexseq_all_pairs_a2.R",
            common_arguments=[
                "--counts", str(counts),
                "--catalog", str(catalog),
                "--samples", str(results / "00_metadata" / "validated_samples.tsv"),
                "--contrasts", str(results / "00_metadata" / "contrasts.tsv"),
                "--outdir", str(outdir),
                "--min-count", "5",
                "--design", str(plan.project["design"]),
                "--fdr", str(reporting["fdr"]),
                "--min-abs-delta-pau", str(reporting["min_abs_delta_pau"]),
            ],
            outdir=outdir,
            log_dir=results / "logs" / "apa_a2" / genome,
            receipt_root=outdir / ".receipts",
            index_path=index_path,
            parallel_jobs=plan.project["resources"]["apa_a2"]["contrast_parallel_jobs"],
            threads=1,
            memory_gb=plan.project["resources"]["apa_a2"]["contrast_memory_gb"],
            output_suffixes=[
                ".apa_a2_sites.tsv",
                ".apa_a2_genes.tsv",
                ".apa_a2_pair_deltas.tsv",
                ".apa_a2_shifts.tsv",
                ".apa_a2_audit.tsv",
            ],
            signature_inputs=[counts, catalog],
            signature_parameters={
                "design": plan.project["design"],
                "genome": genome,
                "reporting": reporting,
                "effect_method": "raw_count_within_gene_pau_v1",
            },
            dry_run=False,
            force=force,
        )
        pcpa_path = output_root / genome / "candidate_pcpa.tsv"
        _write_differential_pcpa(
            active_root / genome / "pcpa_candidate_catalog.tsv",
            index_path,
            pcpa_path,
            float(reporting["fdr"]),
            float(reporting["min_abs_delta_pau"]),
        )
        outputs.extend([index_path, pcpa_path])
        with index_path.open(encoding="utf-8", newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                outputs.extend(
                    Path(row[field])
                    for field in (
                        "result_file",
                        "shift_file",
                        "gene_summary_file",
                        "pair_delta_file",
                        "audit_file",
                    )
                )

    if outputs:
        signature = signature_for(inputs, {
            "module": "apa_a2",
            "contrasts": plan.contrasts,
            "design": plan.project["design"],
            "reporting": reporting,
            "effect_method": "raw_count_within_gene_pau_v1",
        })
        write_receipt(
            "apa_a2_corrected",
            output_root,
            signature,
            outputs,
            ["rna-ends2tracks", "apa-a2"],
        )
    event(
        log_dir,
        "apa_a2",
        "completed",
        "Completed independent DEXSeq tests and corrected raw-count PAU effects",
    )
