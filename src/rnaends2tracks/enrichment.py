from __future__ import annotations

import csv
import math
from pathlib import Path
from typing import Any

from .config import RunPlan, signature_for
from .execution import run_bounded
from .external import event, progress_events, require_tools, run
from .receipts import receipt_valid, write_receipt
from .statistics import r_environment


OUTPUT_NAMES = (
    "ora.tsv", "gsea.tsv", "mapping_audit.tsv", "enrichment.pdf", "enrichment.png", "provenance.tsv",
)
INDEX_FIELDS = (
    "analysis_type", "genome", "contrast_id", "status", "prepared_gene_table",
    "background_genes", "foreground_genes", "ora_terms", "significant_ora_terms",
    "gsea_terms", "significant_gsea_terms", "ora_file", "gsea_file", "plot_pdf",
    "plot_png", "mapping_audit", "provenance_file",
)


def _read(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _write(path: Path, rows: list[dict[str, Any]], fields: tuple[str, ...] | list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fields), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def _number(value: object) -> float | None:
    try:
        result = float(str(value))
    except (TypeError, ValueError):
        return None
    return result if math.isfinite(result) else None


def _rank(row: dict[str, str], effect: float | None, pvalue_column: str = "pvalue") -> float:
    statistic = _number(row.get("stat"))
    if statistic is not None:
        return statistic
    pvalue = _number(row.get(pvalue_column))
    if effect is None or pvalue is None:
        return effect or 0.0
    return math.copysign(-math.log10(max(pvalue, 1e-300)), effect)


def _dge_gene_table(source: Path, target: Path, fdr: float, min_lfc: float) -> tuple[int, int]:
    rows: list[dict[str, Any]] = []
    for row in _read(source):
        gene_id = row.get("gene_id", "")
        if not gene_id:
            continue
        padj = _number(row.get("padj"))
        effect = _number(row.get("log2FoldChange"))
        selected = padj is not None and padj <= fdr and effect is not None and abs(effect) >= min_lfc
        direction = ""
        if selected:
            direction = "any_dge;upregulated" if effect > 0 else "any_dge;downregulated"
        rows.append({
            "gene_id": gene_id, "background": 1, "foreground": int(selected),
            "direction": direction, "rank_score": _rank(row, effect),
        })
    _write(target, rows, ["gene_id", "background", "foreground", "direction", "rank_score"])
    return len(rows), sum(int(row["foreground"]) for row in rows)


def _pcpa_by_gene(path: Path, contrast_id: str) -> dict[str, float]:
    values: dict[str, float] = {}
    if not path.is_file():
        return values
    for row in _read(path):
        if row.get("contrast_id") != contrast_id or not row.get("gene_id"):
            continue
        delta = _number(row.get("delta_PAU"))
        if delta is None:
            continue
        gene = row["gene_id"]
        if gene not in values or abs(delta) > abs(values[gene]):
            values[gene] = delta
    return values


def _apa_gene_table(
    source: Path, pcpa_path: Path, target: Path, contrast_id: str, fdr: float, min_delta: float,
) -> tuple[int, int]:
    pcpa = _pcpa_by_gene(pcpa_path, contrast_id)
    rows: list[dict[str, Any]] = []
    for row in _read(source):
        gene = row.get("gene_id", "")
        if not gene:
            continue
        padj = _number(row.get("gene_padj"))
        delta = _number(row.get("max_abs_delta_PAU"))
        selected_apa = padj is not None and padj <= fdr and delta is not None and delta >= min_delta
        directions: list[str] = []
        if selected_apa:
            directions.append("any_apa")
            shift = row.get("shift", "").lower()
            if shift in {"distal", "proximal"}:
                directions.append(f"{shift}_shift")
        if gene in pcpa:
            directions.extend(["any_pcpa", "pcpa_increased" if pcpa[gene] > 0 else "pcpa_decreased"])
        score = _number(row.get("signed_shift_score"))
        if score is None:
            score = _number(row.get("weighted_transcript_position_shift_nt")) or 0.0
        rows.append({
            "gene_id": gene, "background": 1, "foreground": int(bool(directions)),
            "direction": ";".join(dict.fromkeys(directions)), "rank_score": score,
        })
    _write(target, rows, ["gene_id", "background", "foreground", "direction", "rank_score"])
    return len(rows), sum(int(row["foreground"]) for row in rows)


def _count_significant(path: Path, fdr: float) -> tuple[int, int]:
    rows = _read(path)
    return len(rows), sum((_number(row.get("padj")) or math.inf) <= fdr for row in rows)


def _index_row(
    analysis_type: str, genome: str, contrast_id: str, outdir: Path,
    gene_table: Path, background: int, foreground: int, fdr: float,
) -> dict[str, Any]:
    ora_total, ora_significant = _count_significant(outdir / "ora.tsv", fdr)
    gsea_total, gsea_significant = _count_significant(outdir / "gsea.tsv", fdr)
    return {
        "analysis_type": analysis_type, "genome": genome, "contrast_id": contrast_id,
        "status": "PASS", "prepared_gene_table": str(gene_table),
        "background_genes": background, "foreground_genes": foreground,
        "ora_terms": ora_total, "significant_ora_terms": ora_significant,
        "gsea_terms": gsea_total, "significant_gsea_terms": gsea_significant,
        "ora_file": str(outdir / "ora.tsv"), "gsea_file": str(outdir / "gsea.tsv"),
        "plot_pdf": str(outdir / "enrichment.pdf"), "plot_png": str(outdir / "enrichment.png"),
        "mapping_audit": str(outdir / "mapping_audit.tsv"),
        "provenance_file": str(outdir / "provenance.tsv"),
    }


def _ensure_apa_a_gene_summaries(
    plan: RunPlan, results: Path, script_root: Path,
) -> None:
    """Upgrade pre-alpha.10 APA-A indexes before enrichment.

    Older result indexes did not publish the gene-level summary required by
    enrichment.  Re-running only the APA-A statistical layer is exact and much
    safer than approximating gene-level q-values from site-level tables.
    """
    legacy: list[Path] = []
    for genome in plan.references:
        index = results / "06_apa_a_mcell2019" / genome / "dexseq" / "result_index.tsv"
        rows = _read(index) if index.is_file() else []
        if rows and any(
            not row.get("gene_summary_file")
            or not Path(row.get("gene_summary_file", "")).is_file()
            for row in rows
        ):
            legacy.append(index)
    if not legacy:
        return
    event(
        results / "logs", "enrichment", "progress",
        "Detected pre-alpha.10 APA-A indexes; regenerating the APA-A statistical layer "
        "to publish exact gene-level enrichment sources",
    )
    from .apa_mcell import apa_statistics_stage

    apa_statistics_stage(plan, results, script_root, dry_run=False, force=True)
    unresolved = [
        path for path in legacy
        if any(
            not row.get("gene_summary_file")
            or not Path(row.get("gene_summary_file", "")).is_file()
            for row in _read(path)
        )
    ]
    if unresolved:
        raise RuntimeError(
            "APA-A compatibility regeneration did not publish gene summaries: "
            + ", ".join(map(str, unresolved))
        )


def enrichment(
    plan: RunPlan, results: Path, script_root: Path, dry_run: bool = False, force: bool = False,
) -> None:
    modules = plan.project.get("modules", {})
    dge_enabled = bool(modules.get("dge_enrichment", True) and modules.get("gene_expression", True))
    apa_enabled = bool(modules.get("apa_enrichment", True) and modules.get("apa_a", True))
    apa_b_enabled = bool(modules.get("apa_enrichment", True) and plan.project.get("apa_b", {}).get("enabled", False))
    log_dir = results / "logs"
    outroot = results / "10_reports" / "enrichment_summary"
    if not (dge_enabled or apa_enabled or apa_b_enabled):
        event(log_dir, "enrichment", "disabled", "DGE and APA enrichment are disabled")
        return
    if dry_run:
        event(log_dir, "enrichment", "dry_run", "Would run ORA and ranked GSEA for enabled DGE/APA contrasts")
        return
    require_tools(["Rscript"])
    if apa_enabled:
        _ensure_apa_a_gene_summaries(plan, results, script_root)
    settings = plan.project["enrichment"]
    fdr = float(plan.project["reporting"]["fdr"])
    enrichment_padj = float(settings["padj"])
    jobs: list[tuple[str, Any]] = []
    signature_inputs: list[Path] = []

    def add_job(analysis_type: str, genome: str, species: str, contrast_id: str, source: Path, pcpa: Path | None) -> None:
        if not source.is_file():
            raise RuntimeError(f"Missing {analysis_type} enrichment source: {source}")
        if analysis_type == "dge":
            job_dir = results / "05_gene_expression" / genome / "enrichment" / contrast_id
        elif analysis_type == "apa_a":
            job_dir = results / "06_apa_a_mcell2019" / genome / "enrichment" / contrast_id
        else:
            job_dir = results / "07_apa_b" / genome / "enrichment" / contrast_id
        gene_table = job_dir / "prepared_gene_table.tsv"
        inputs = [source] + ([pcpa] if pcpa is not None and pcpa.is_file() else [])
        signature_inputs.extend(inputs)
        job_signature = signature_for(inputs, {
            "module": "enrichment", "analysis_type": analysis_type, "genome": genome,
            "contrast_id": contrast_id, "settings": settings, "fdr": fdr,
        })

        def worker() -> dict[str, Any]:
            if not force and receipt_valid(job_dir, job_signature):
                rows = _read(job_dir / "job_index.tsv")
                if rows:
                    return rows[0]
            job_dir.mkdir(parents=True, exist_ok=True)
            if analysis_type == "dge":
                background, foreground = _dge_gene_table(
                    source, gene_table, fdr, float(settings["dge_min_abs_lfc"]),
                )
            else:
                background, foreground = _apa_gene_table(
                    source, pcpa or Path(""), gene_table, contrast_id, fdr,
                    float(settings["apa_min_abs_delta_pau"]),
                )
            command = [
                "Rscript", str(script_root / "R" / "gene_set_enrichment.R"),
                "--input", str(gene_table), "--outdir", str(job_dir),
                "--analysis-type", analysis_type, "--species", species,
                "--genome", genome, "--contrast-id", contrast_id,
                "--ora", str(settings["ora"]).lower(), "--gsea", str(settings["gsea"]).lower(),
                "--go", str(settings["go"]).lower(), "--reactome", str(settings["reactome"]).lower(),
                "--hallmarks", str(settings["hallmarks"]).lower(),
                "--padj", str(enrichment_padj), "--min-size", str(settings["min_geneset_size"]),
                "--max-size", str(settings["max_geneset_size"]),
            ]
            resource = plan.project["resources"]["enrichment"]
            run(command, log_dir / "enrichment" / genome / contrast_id / f"{analysis_type}.log",
                env=r_environment(resource["threads"], resource["memory_gb"]))
            row = _index_row(
                analysis_type, genome, contrast_id, job_dir, gene_table, background, foreground, enrichment_padj,
            )
            job_index = job_dir / "job_index.tsv"
            _write(job_index, [row], INDEX_FIELDS)
            write_receipt("enrichment_job", job_dir, job_signature,
                          [gene_table, *(job_dir / name for name in OUTPUT_NAMES), job_index],
                          ["rna-ends2tracks", "enrichment", analysis_type, contrast_id])
            return row

        jobs.append((f"{analysis_type}:{genome}:{contrast_id}", worker))

    for genome, reference in plan.references.items():
        species = str(reference["species"])
        if dge_enabled:
            for row in _read(results / "05_gene_expression" / genome / "C4_primary_deseq2" / "result_index.tsv"):
                add_job("dge", genome, species, row["contrast_id"], Path(row["result_file"]), None)
        if apa_enabled:
            pcpa = results / "06_apa_a_mcell2019" / genome / "candidate_pcpa.tsv"
            for row in _read(results / "06_apa_a_mcell2019" / genome / "dexseq" / "result_index.tsv"):
                add_job("apa_a", genome, species, row["contrast_id"], Path(row["gene_summary_file"]), pcpa)
        if apa_b_enabled:
            pcpa = results / "07_apa_b" / genome / "candidate_pcpa.tsv"
            for row in _read(results / "07_apa_b" / genome / "drimseq" / "result_index.tsv"):
                add_job("apa_b", genome, species, row["contrast_id"], Path(row["gene_summary_file"]), pcpa)
    if not jobs:
        raise RuntimeError("Enrichment is enabled, but no upstream contrast outputs were found")
    resource = plan.project["resources"]["enrichment"]
    rows = run_bounded(
        "enrichment", jobs, resource["parallel_jobs"], results / ".checkpoints" / "timings" / "enrichment",
        progress=progress_events(log_dir, "enrichment", len(jobs), "analysis"),
    )
    index = outroot / "enrichment_index.tsv"
    _write(index, rows, INDEX_FIELDS)
    module_signature = signature_for(signature_inputs, {
        "module": "enrichment", "settings": settings,
        "dge_enabled": dge_enabled, "apa_enabled": apa_enabled, "apa_b_enabled": apa_b_enabled,
    })
    outputs = [index]
    for row in rows:
        outputs.extend([Path(row[field]) for field in (
            "prepared_gene_table", "ora_file", "gsea_file", "plot_pdf", "plot_png", "mapping_audit", "provenance_file",
        )])
    write_receipt("enrichment", outroot, module_signature, outputs, ["rna-ends2tracks", "enrichment"])
    event(log_dir, "enrichment", "completed", f"Completed {len(rows)} DGE/APA enrichment analyses")
