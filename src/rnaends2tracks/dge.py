from __future__ import annotations

import csv
import math
from pathlib import Path

from .config import RunPlan, signature_for
from .external import event, require_tools, run
from .receipts import receipt_valid, write_receipt
from .statistics import r_environment, run_r_contrasts


def _read_matrix(path: Path, id_column: str) -> tuple[list[str], dict[str, list[int]]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader((line for line in handle if not line.startswith("#")), delimiter="\t")
        if id_column not in (reader.fieldnames or []):
            raise RuntimeError(f"Missing {id_column} in {path}")
        columns = [column for column in (reader.fieldnames or []) if column != id_column]
        values = {row[id_column]: [int(float(row[column])) for column in columns] for row in reader}
    return columns, values


def _correlation(left: list[int], right: list[int]) -> float:
    x = [math.log2(value + 1) for value in left]
    y = [math.log2(value + 1) for value in right]
    if len(x) < 2:
        return float("nan")
    x_mean, y_mean = sum(x) / len(x), sum(y) / len(y)
    numerator = sum((a - x_mean) * (b - y_mean) for a, b in zip(x, y))
    denominator = math.sqrt(sum((a - x_mean) ** 2 for a in x) * sum((b - y_mean) ** 2 for b in y))
    return numerator / denominator if denominator else float("nan")


def _write_c4_c5_diagnostic(
    c4: Path, c5: Path, sample_ids: list[str], output: Path, discrepancy_output: Path,
) -> None:
    _, c4_rows = _read_matrix(c4, "gene_id")
    with c5.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader((line for line in handle if not line.startswith("#")), delimiter="\t")
        annotation = {"Geneid", "Chr", "Start", "End", "Strand", "Length"}
        count_columns = [column for column in (reader.fieldnames or []) if column not in annotation]
        c5_rows = {row["Geneid"]: [int(float(row[column])) for column in count_columns] for row in reader}
    shared = sorted(set(c4_rows) & set(c5_rows))
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["sample_id", "shared_genes", "pearson_log2p1", "C4_total", "C5_total"])
        for index, sample_id in enumerate(sample_ids):
            left = [c4_rows[gene][index] for gene in shared]
            right = [c5_rows[gene][index] for gene in shared]
            writer.writerow([sample_id, len(shared), _correlation(left, right), sum(left), sum(right)])
    with discrepancy_output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["sample_id", "gene_id", "C4_count", "C5_count", "C4_CPM", "C5_CPM", "log2_CPM_ratio"])
        for index, sample_id in enumerate(sample_ids):
            c4_total = sum(values[index] for values in c4_rows.values())
            c5_total = sum(values[index] for values in c5_rows.values())
            for gene in shared:
                left, right = c4_rows[gene][index], c5_rows[gene][index]
                left_cpm = 1_000_000 * left / c4_total if c4_total else 0.0
                right_cpm = 1_000_000 * right / c5_total if c5_total else 0.0
                ratio = math.log2((left_cpm + 0.5) / (right_cpm + 0.5))
                if abs(ratio) >= 2:
                    writer.writerow([sample_id, gene, left, right, left_cpm, right_cpm, ratio])


def gene_expression(plan: RunPlan, results: Path, script_root: Path, dry_run: bool = False, force: bool = False) -> None:
    module_dir = results / "05_gene_expression"
    log_dir = results / "logs"
    module_dir.mkdir(parents=True, exist_ok=True)
    if not plan.project.get("modules", {}).get("gene_expression", True):
        event(log_dir, "gene_expression", "disabled", "RUN_GENE_EXPRESSION=false")
        return
    signature_inputs = [
        results / "02_alignment" / sample["sample_id"] / f"{sample['sample_id']}.bam" for sample in plan.samples
    ] + [
        results / "04_active_pas" / genome / "C4_active_pas_gene_counts.tsv"
        for genome in (plan.references or {plan.reference["assembly"]: plan.reference})
        if any(sample["genome"] == genome for sample in plan.samples)
    ] + [
        Path(plan.reference_for(genome)["gtf"])
        for genome in (plan.references or {plan.reference["assembly"]: plan.reference})
        if any(sample["genome"] == genome for sample in plan.samples)
    ]
    signature = signature_for(signature_inputs, {
        "module": "gene_expression", "contrasts": plan.contrasts,
        "design": plan.project["design"], "reporting": plan.project["reporting"],
    }) if not dry_run else "dry-run"
    if not force and not dry_run and receipt_valid(module_dir, signature):
        event(log_dir, "gene_expression", "skipped", "Valid C4/C5 gene-expression receipt")
        return
    if not dry_run:
        require_tools(["featureCounts", "Rscript"])
    expected: list[Path] = []
    for genome in plan.references or {plan.reference["assembly"]: plan.reference}:
        samples = [sample for sample in plan.samples if sample["genome"] == genome]
        contrasts = [contrast for contrast in plan.contrasts if contrast.get("genome", genome) == genome]
        if not samples:
            continue
        reference = plan.reference_for(genome)
        sample_ids = [sample["sample_id"] for sample in samples]
        bams = [results / "02_alignment" / sample_id / f"{sample_id}.bam" for sample_id in sample_ids]
        c4 = results / "04_active_pas" / genome / "C4_active_pas_gene_counts.tsv"
        outdir = module_dir / genome
        c5 = outdir / "C5_featureCounts_diagnostic.tsv"
        diagnostic = outdir / "C4_vs_C5_library_correlations.tsv"
        discrepancies = outdir / "C4_vs_C5_large_discrepancies.tsv"
        primary = outdir / "C4_primary_deseq2"
        factor_path = results / "04_active_pas" / genome / "C4_track_size_factors.tsv"
        index_path = primary / "result_index.tsv"
        resource = plan.project["resources"]["dge"]
        run([
            "featureCounts", "-T", str(resource["featurecounts_threads"]), "-a", reference["gtf"],
            "-o", str(c5), "-t", "exon", "-g", "gene_id", "-s", "2", "--primary", *map(str, bams),
        ], log_dir / "gene_expression" / genome / "featureCounts.log", dry_run)
        common = [
            "--counts", str(c4), "--samples", str(results / "00_metadata" / "validated_samples.tsv"),
            "--contrasts", str(results / "00_metadata" / "contrasts.tsv"), "--design", str(plan.project["design"]),
            "--outdir", str(primary), "--fdr", str(plan.project["reporting"]["fdr"]),
            "--factor-output", str(factor_path),
        ]
        run([
            "Rscript", str(script_root / "R" / "deseq2_c4.R"), *common, "--mode", "qc",
        ], log_dir / "gene_expression" / genome / "C4_qc.log", dry_run, env=r_environment(1, resource["contrast_memory_gb"]))
        if not dry_run:
            factor_receipt_dir = factor_path.parent / ".track_factor_receipt"
            factor_receipt_dir.mkdir(parents=True, exist_ok=True)
            factor_signature = signature_for(
                [c4, results / "00_metadata" / "validated_samples.tsv"],
                {"module": "C4_track_size_factors", "genome": genome, "estimator": "DESeq2_poscounts"},
            )
            write_receipt(
                "C4_track_size_factors", factor_receipt_dir,
                factor_signature, [factor_path], ["rna-ends2tracks", "C4_track_size_factors", genome],
            )
        if contrasts:
            genome_plan = RunPlan(plan.project, samples, [row for row in plan.sample_rows if row["genome"] == genome], contrasts, reference, {genome: reference})
            run_r_contrasts(
                module=f"dge_{genome}", plan=genome_plan, results=results,
                script=script_root / "R" / "deseq2_c4.R", common_arguments=[*common, "--mode", "contrast"],
                outdir=primary, log_dir=log_dir / "gene_expression" / genome / "contrasts",
                receipt_root=primary / ".receipts", index_path=index_path,
                parallel_jobs=resource["contrast_parallel_jobs"], threads=resource["contrast_threads"],
                memory_gb=resource["contrast_memory_gb"],
                output_suffixes=[".deseq2.tsv", ".deseq2_model.rds", ".MA.pdf"],
                signature_inputs=[c4, results / "00_metadata" / "validated_samples.tsv", results / "00_metadata" / "contrasts.tsv"],
                signature_parameters={
                    "genome": genome, "design": plan.project["design"],
                    "reporting": plan.project["reporting"],
                }, dry_run=dry_run, force=force,
            )
        if not dry_run:
            _write_c4_c5_diagnostic(c4, c5, sample_ids, diagnostic, discrepancies)
        expected.extend([c5, Path(str(c5) + ".summary"), diagnostic, discrepancies,
                         factor_path, primary / "C4_deseq2_model.rds",
                         primary / "C4_normalized_counts.tsv", primary / "C4_vst_pca.pdf"])
        if contrasts:
            expected.append(index_path)
            if not dry_run:
                with index_path.open(encoding="utf-8", newline="") as handle:
                    for row in csv.DictReader(handle, delimiter="\t"):
                        contrast_id = row["contrast_id"]
                        expected.extend([Path(row["result_file"]),
                                         primary / f"{contrast_id}.deseq2_model.rds",
                                         primary / f"{contrast_id}.MA.pdf"])
    if not dry_run:
        write_receipt("gene_expression", module_dir, signature, expected, ["rna-ends2tracks", "gene_expression"])
    event(log_dir, "gene_expression", "dry_run" if dry_run else "completed", "C4 primary DGE; C5 diagnostic only")
