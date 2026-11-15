"""PAS-independent gene counting and paired effect-hypothesis analyses."""
from __future__ import annotations

import csv
from pathlib import Path

from .config import RunPlan, signature_for
from .expression_methods import NEW_METHOD, sample_bam
from .external import event, require_tools, run
from .gex_annotation import build_annotation
from .receipts import receipt_valid, sha256, write_receipt
from .statistics import r_environment, run_r_contrasts


def export_counts(source: Path, target: Path, sample_ids: list[str], bams: list[Path]) -> None:
    """Validate featureCounts columns explicitly; never silently reorder samples."""
    with source.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader((line for line in handle if not line.startswith("#")), delimiter="\t")
        columns = list(reader.fieldnames or [])
        if columns[:6] != ["Geneid", "Chr", "Start", "End", "Strand", "Length"]:
            raise RuntimeError("Unexpected featureCounts annotation columns")
        counts = columns[6:]
        if len(counts) != len(bams) or any(Path(column).resolve() != bam.resolve() for column, bam in zip(counts, bams)):
            raise RuntimeError("featureCounts BAM columns do not match the requested samples/order")
        seen = set()
        with target.open("w", encoding="utf-8", newline="") as output:
            writer = csv.writer(output, delimiter="\t", lineterminator="\n")
            writer.writerow(["gene_id", *sample_ids])
            for row in reader:
                gene = row["Geneid"]
                if not gene or gene in seen:
                    raise RuntimeError("Empty or duplicate featureCounts gene ID")
                seen.add(gene)
                values = [int(row[column]) for column in counts]
                if any(value < 0 for value in values):
                    raise RuntimeError("Negative gene count")
                writer.writerow([gene, *values])
        if not seen:
            raise RuntimeError("featureCounts produced no genes")


def extended_gene_expression(plan: RunPlan, results: Path, script_root: Path, dry_run=False, force=False) -> list[Path]:
    settings = plan.project["gene_expression"]
    resource = plan.project["resources"]["dge"]
    script = script_root / "R" / "deseq2_gex.R"
    metadata = [results / "00_metadata" / name for name in ("validated_samples.tsv", "contrasts.tsv")]
    expected = []
    if not dry_run:
        require_tools(["featureCounts", "Rscript"])
    for genome, reference in (plan.references or {plan.reference["assembly"]: plan.reference}).items():
        samples = [row for row in plan.samples if row["genome"] == genome]
        if not samples:
            continue
        contrasts = [row for row in plan.contrasts if row.get("genome", genome) == genome]
        outdir = results / "05_gene_expression" / genome / NEW_METHOD
        outdir.mkdir(parents=True, exist_ok=True)
        temporary = Path(plan.project.get("tmp_dir") or outdir / ".tmp") / "gene_expression" / genome
        temporary.mkdir(parents=True, exist_ok=True)
        tool_env = {"TMPDIR": str(temporary), "TMP": str(temporary), "TEMP": str(temporary)}
        logs = results / "logs" / "gene_expression" / genome / NEW_METHOD
        annotation_dir = outdir / "annotation"
        annotation = {"saf": annotation_dir / "extended_exons.saf"}
        if not dry_run:
            annotation = build_annotation(Path(reference["gtf"]), Path(reference["chrom_sizes"]), annotation_dir,
                                          settings["extension_max_nt"], settings["tss_buffer_nt"])
        bams = [sample_bam(plan.project, results, row["sample_id"]) for row in samples]
        featurecounts = outdir / "featureCounts.tsv"
        counts = outdir / "raw_counts.tsv"
        count_outputs = [featurecounts, Path(str(featurecounts)+".summary"), counts]
        paired = ["-p", "--countReadPairs"] if samples[0].get("library_layout") == "PE" else []
        command = ["featureCounts", "-T", str(resource["featurecounts_threads"]), "-F", "SAF", "-a", str(annotation["saf"]),
                   "-o", str(featurecounts), "-s", "2", "--primary", "--tmpDir", str(temporary), *paired, *map(str, bams)]
        count_signature = "dry-run" if dry_run else signature_for([*bams, annotation["manifest"]],
            {"command": command, "samples": [row["sample_id"] for row in samples], "counter_schema": 1})
        if force or dry_run or not receipt_valid(outdir / ".count_receipt", count_signature):
            run(command, logs / "featureCounts.log", dry_run, env=tool_env)
            if not dry_run:
                export_counts(featurecounts, counts, [row["sample_id"] for row in samples], bams)
                write_receipt("extended_exon_counts", outdir / ".count_receipt", count_signature, count_outputs, command)
        else:
            event(results / "logs", "gene_expression", "skipped", f"Reusing {genome} extended-exon raw counts")
        common = ["--counts", str(counts), "--samples", str(metadata[0]), "--contrasts", str(metadata[1]),
                  "--design", str(plan.project["design"]), "--fdr", str(plan.project["reporting"]["fdr"]),
                  "--min-count", str(settings["min_count"]), "--min-total-count", str(settings["min_total_count"]),
                  "--lfc-threshold", str(settings["lfc_threshold"]), "--stringent", str(settings["stringent"]).lower(),
                  "--relaxed", str(settings["relaxed"]).lower(), "--model-code-hash", sha256(script), "--force-fit", str(force).lower()]
        qc_outputs = [outdir / "qc" / name for name in (
            "normalized_counts.tsv", "library_sizes.tsv", "pca.tsv", "pca.png", "pca.pdf",
            "sample_distances.tsv", "sample_distances.png", "sample_distances.pdf", "sessionInfo.txt")]
        qc_signature = "dry-run" if dry_run else signature_for([counts, *metadata, script],
            {"design": plan.project["design"], "min_count": settings["min_count"], "min_total": settings["min_total_count"]})
        if force or dry_run or not receipt_valid(outdir / "qc", qc_signature):
            run(["Rscript", str(script), *common, "--mode", "qc", "--outdir", str(outdir)], logs / "qc.log", dry_run,
                env={**r_environment(1, resource["contrast_memory_gb"]), **tool_env})
            if not dry_run:
                write_receipt("extended_exon_qc", outdir / "qc", qc_signature, qc_outputs, ["rna-ends2tracks", "gex_qc"])
        contrast_dir = outdir / "contrasts"
        suffixes = ["/model.rds", "/eligibility.tsv", "/filter_parameters.tsv", "/test_universe.tsv", "/ranking.tsv",
                    "/normalized_counts.tsv", "/sessionInfo.txt"]
        for mode, enabled in (("stringent", settings["stringent"]), ("relaxed_fdr", settings["relaxed"])):
            if enabled:
                suffixes.extend(f"/{mode}/{name}" for name in ("deseq2.tsv", "significant_genes.tsv", "MA.png", "MA.pdf", "volcano.png", "volcano.pdf"))
        index = outdir / "result_index.tsv"
        if contrasts:
            genome_plan = RunPlan(plan.project, samples, [row for row in plan.sample_rows if row["genome"] == genome], contrasts, reference, {genome: reference})
            run_r_contrasts(module=f"gex_{genome}", plan=genome_plan, results=results, script=script,
                common_arguments=[*common, "--mode", "contrast", "--outdir", str(contrast_dir)], outdir=contrast_dir,
                log_dir=logs / "contrasts", receipt_root=outdir / ".receipts", index_path=index,
                parallel_jobs=resource["contrast_parallel_jobs"], threads=resource["contrast_threads"], memory_gb=resource["contrast_memory_gb"],
                output_suffixes=suffixes, signature_inputs=[counts, *metadata, script],
                signature_parameters={"settings": settings, "fdr": plan.project["reporting"]["fdr"]}, dry_run=dry_run, force=force, environment=tool_env)
        elif not dry_run:
            index.write_text("expression_method\tanalysis_mode\tcontrast_id\tresult_file\n", encoding="utf-8")
        outputs = [*annotation.values(), *count_outputs, *qc_outputs, index]
        outputs += [contrast_dir / (row["contrast_id"]+suffix) for row in contrasts for suffix in suffixes]
        expected.extend(outputs)
        if not dry_run:
            signature = signature_for([*metadata, script, annotation["manifest"]], {"settings": settings, "counts": count_signature})
            write_receipt("extended_gene_expression", outdir, signature, outputs, ["rna-ends2tracks", NEW_METHOD])
    event(results / "logs", "gene_expression", "dry_run" if dry_run else "completed", "Extended-exon counts; stringent and relaxed DE; shared ordinary-Wald GSEA ranking")
    return expected
