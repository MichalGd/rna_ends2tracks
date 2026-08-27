from __future__ import annotations

import csv
import html
import json
from collections import Counter
from pathlib import Path
from typing import Any

from .config import RunPlan, signature_for, workflow_requirements
from .external import event
from .receipts import receipt_valid, write_receipt


CONTRAST_SUMMARY_FIELDS = [
    "genome", "contrast_id", "numerator", "denominator", "design_mode",
    "dge_tested_genes", "dge_significant", "dge_up", "dge_down",
    "apa_a_tested_sites", "apa_a_significant_sites", "apa_a_distal_genes",
    "apa_a_proximal_genes", "apa_a_pcpa", "apa_b_tested_sites",
    "apa_b_confirmed_sites", "apa_b_pcpa", "apa_compared_sites",
    "apa_direction_agree", "apa_direction_disagree",
]


def _rows(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        return []
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _receipt_status(path: Path) -> str:
    try:
        payload = json.loads((path / "run_receipt.json").read_text(encoding="utf-8"))
        return "PASS" if payload.get("exit_status") == 0 else "FAIL"
    except (OSError, ValueError):
        return "NOT_RUN"


def _number(value: str | None) -> float | None:
    if value in {None, "", "NA", "NaN", "nan"}:
        return None
    try:
        return float(value)
    except ValueError:
        return None


def _index_by_contrast(path: Path) -> dict[str, dict[str, str]]:
    return {row["contrast_id"]: row for row in _rows(path) if row.get("contrast_id")}


def _count_pcpa(path: Path) -> Counter[str]:
    return Counter(row.get("contrast_id", "") for row in _rows(path) if row.get("contrast_id"))


def _required_rows(path: Path, label: str) -> list[dict[str, str]]:
    if not path.is_file():
        raise RuntimeError(f"Missing {label} report source: {path}")
    return _rows(path)


def _dge_counts(path: Path, fdr: float) -> dict[str, int]:
    rows = _required_rows(path, "DGE")
    significant = [
        row for row in rows
        if (_number(row.get("padj")) is not None and _number(row.get("padj")) <= fdr)
    ]
    return {
        "dge_tested_genes": len(rows),
        "dge_significant": len(significant),
        "dge_up": sum((_number(row.get("log2FoldChange")) or 0) > 0 for row in significant),
        "dge_down": sum((_number(row.get("log2FoldChange")) or 0) < 0 for row in significant),
    }


def _tested_significant(path: Path, adjusted_column: str, fdr: float) -> tuple[int, int]:
    rows = _required_rows(path, adjusted_column)
    return len(rows), sum(
        _number(row.get(adjusted_column)) is not None
        and _number(row.get(adjusted_column)) <= fdr
        for row in rows
    )


def _require_index_count(
    module: str, contrast_id: str, index: dict[str, str], field: str, observed: int,
) -> None:
    declared = index.get(field, "")
    if declared and int(float(declared)) != observed:
        raise RuntimeError(
            f"{module} summary mismatch for {contrast_id}: index {field}={declared}, observed={observed}"
        )


def _comparison_counts(path: Path) -> dict[str, dict[str, int]]:
    counts: dict[str, dict[str, int]] = {}
    for row in _rows(path):
        contrast = row.get("contrast_id", "")
        if not contrast:
            continue
        record = counts.setdefault(contrast, {
            "apa_compared_sites": 0, "apa_direction_agree": 0, "apa_direction_disagree": 0,
        })
        record["apa_compared_sites"] += 1
        agreement = str(row.get("direction_agreement", "")).lower()
        if agreement == "true":
            record["apa_direction_agree"] += 1
        elif agreement == "false":
            record["apa_direction_disagree"] += 1
    return counts


def _contrast_summary(
    plan: RunPlan, results: Path, contrasts: list[dict[str, str]],
) -> list[dict[str, Any]]:
    """Create one independently testable row per configured contrast."""
    fdr = float(plan.project.get("reporting", {}).get("fdr", 0.05))
    comparison = _comparison_counts(results / "08_apa_comparison" / "effect_concordance.tsv")
    summary: list[dict[str, Any]] = []
    for genome in sorted({row.get("genome", "") for row in contrasts}):
        genome_contrasts = [row for row in contrasts if row.get("genome", "") == genome]
        dge_index = _index_by_contrast(
            results / "05_gene_expression" / genome / "C4_primary_deseq2" / "result_index.tsv"
        )
        apa_a_index = _index_by_contrast(
            results / "06_apa_a_mcell2019" / genome / "dexseq" / "result_index.tsv"
        )
        apa_b_index = _index_by_contrast(
            results / "07_apa_b" / genome / "drimseq" / "result_index.tsv"
        )
        apa_a_pcpa_path = results / "06_apa_a_mcell2019" / genome / "candidate_pcpa.tsv"
        apa_b_pcpa_path = results / "07_apa_b" / genome / "candidate_pcpa.tsv"
        apa_a_pcpa = _count_pcpa(apa_a_pcpa_path) if apa_a_pcpa_path.is_file() else None
        apa_b_pcpa = _count_pcpa(apa_b_pcpa_path) if apa_b_pcpa_path.is_file() else None
        for contrast in genome_contrasts:
            contrast_id = contrast["contrast_id"]
            row: dict[str, Any] = {
                "genome": genome,
                "contrast_id": contrast_id,
                "numerator": contrast.get("numerator", ""),
                "denominator": contrast.get("denominator", ""),
                "design_mode": contrast.get("design_mode", ""),
                "dge_tested_genes": "",
                "dge_significant": "",
                "dge_up": "",
                "dge_down": "",
                "apa_a_tested_sites": "",
                "apa_a_significant_sites": "",
                "apa_a_distal_genes": "",
                "apa_a_proximal_genes": "",
                "apa_a_pcpa": "" if apa_a_pcpa is None else apa_a_pcpa[contrast_id],
                "apa_b_tested_sites": "",
                "apa_b_confirmed_sites": "",
                "apa_b_pcpa": "" if apa_b_pcpa is None else apa_b_pcpa[contrast_id],
                "apa_compared_sites": "",
                "apa_direction_agree": "",
                "apa_direction_disagree": "",
            }
            if contrast_id in dge_index:
                index = dge_index[contrast_id]
                counts = _dge_counts(Path(index["result_file"]), fdr)
                _require_index_count(
                    "DGE", contrast_id, index, "significant", counts["dge_significant"],
                )
                row.update(counts)
            if contrast_id in apa_a_index:
                index = apa_a_index[contrast_id]
                tested, significant = _tested_significant(Path(index["result_file"]), "padj", fdr)
                _require_index_count("APA-A", contrast_id, index, "tested_sites", tested)
                _require_index_count("APA-A", contrast_id, index, "significant_sites", significant)
                shifts = Counter(
                    item.get("shift", "")
                    for item in _required_rows(Path(index["shift_file"]), "APA-A shift")
                )
                row.update({
                    "apa_a_tested_sites": tested,
                    "apa_a_significant_sites": significant,
                    "apa_a_distal_genes": shifts["distal"],
                    "apa_a_proximal_genes": shifts["proximal"],
                })
            if contrast_id in apa_b_index:
                index = apa_b_index[contrast_id]
                tested, confirmed = _tested_significant(
                    Path(index["result_file"]), "stageR_adjusted", fdr,
                )
                _require_index_count("APA-B", contrast_id, index, "tested_sites", tested)
                _require_index_count("APA-B", contrast_id, index, "confirmed_sites", confirmed)
                row.update({
                    "apa_b_tested_sites": tested,
                    "apa_b_confirmed_sites": confirmed,
                })
            row.update(comparison.get(contrast_id, {}))
            summary.append(row)
    return summary


def _write_tsv(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def _html_table(rows: list[dict[str, Any]], fields: list[str]) -> str:
    headings = "".join(f"<th>{html.escape(field.replace('_', ' '))}</th>" for field in fields)
    body = "".join(
        "<tr>" + "".join(f"<td>{html.escape(str(row.get(field, '')))}</td>" for field in fields) + "</tr>"
        for row in rows
    )
    return f"<div class='table-wrap'><table><thead><tr>{headings}</tr></thead><tbody>{body}</tbody></table></div>"


def _browser_assets(results: Path, outdir: Path) -> list[Path]:
    bigwigs = sorted((results / "09_tracks").rglob("*.bw")) if (results / "09_tracks").is_dir() else []
    trackdb = outdir / "UCSC_trackDb.txt"
    with trackdb.open("w", encoding="utf-8") as handle:
        for index, path in enumerate(bigwigs, start=1):
            name = f"rna_ends_{index:04d}"
            relative = path.relative_to(results).as_posix()
            handle.write(f"track {name}\ntype bigWig\nbigDataUrl ../{relative}\nshortLabel {path.stem[:17]}\n")
            handle.write(f"longLabel {relative}\nvisibility full\n\n")
    igv = outdir / "IGV_session.xml"
    resources = "\n".join(
        f'    <Resource path="../{html.escape(path.relative_to(results).as_posix())}"/>' for path in bigwigs
    )
    igv.write_text(f'<?xml version="1.0" encoding="UTF-8"?>\n<Session version="8">\n  <Resources>\n{resources}\n  </Resources>\n</Session>\n',
                   encoding="utf-8")
    return [trackdb, igv]


def make_report(
    plan: RunPlan, results: Path, dry_run: bool = False, force: bool = False,
) -> Path:
    outdir = results / "10_reports"
    log_dir = results / "logs"
    if dry_run:
        event(log_dir, "report", "dry_run", "Would generate Markdown/HTML/TSV and browser-session assets")
        return outdir / "report.html"
    outdir.mkdir(parents=True, exist_ok=True)
    inputs = [results / "00_metadata" / name for name in (
        "contrasts.tsv", "validated_samples.tsv", "resolved_config.json", "warnings.tsv")]
    inputs.extend(
        results / directory / "run_receipt.json"
        for directory in (
            "02_alignment", "03_exact_ends", "04_active_pas", "05_gene_expression",
            "06_apa_a_mcell2019", "07_apa_b", "08_apa_comparison", "09_tracks",
        )
    )
    result_indexes = sorted((results / "05_gene_expression").rglob("result_index.tsv"))
    result_indexes.extend(sorted((results / "06_apa_a_mcell2019").rglob("result_index.tsv")))
    result_indexes.extend(sorted((results / "07_apa_b").rglob("result_index.tsv")))
    inputs.extend(result_indexes)
    for index_path in result_indexes:
        for row in _rows(index_path):
            for field in ("result_file", "shift_file"):
                if row.get(field):
                    inputs.append(Path(row[field]))
    inputs.extend(sorted((results / "06_apa_a_mcell2019").rglob("candidate_pcpa.tsv")))
    inputs.extend(sorted((results / "07_apa_b").rglob("candidate_pcpa.tsv")))
    inputs.append(results / "08_apa_comparison" / "effect_concordance.tsv")
    inputs = [path for path in inputs if path.is_file()]
    signature = signature_for(inputs, {"module": "report", "project": plan.project["project_id"]})
    html_target = outdir / "report.html"
    if not force and receipt_valid(outdir, signature):
        event(log_dir, "report", "skipped", "Valid matching receipt")
        return html_target

    contrasts = _rows(results / "00_metadata" / "contrasts.tsv")
    warnings = _rows(results / "00_metadata" / "warnings.tsv")
    enabled_modules = plan.project.get("modules", {})
    requirements = workflow_requirements(plan.project)
    modules = [
        ("alignment", "02_alignment", True),
        ("exact ends", "03_exact_ends", requirements["exact_ends"]),
        ("active PAS", "04_active_pas", requirements["active_pas"]),
        ("gene expression", "05_gene_expression", enabled_modules.get("gene_expression", True)),
        ("APA-A", "06_apa_a_mcell2019", enabled_modules.get("apa_a", True)),
        ("APA-B", "07_apa_b", plan.project.get("apa_b", {}).get("enabled", False)),
        ("APA comparison", "08_apa_comparison", requirements["apa_comparison"]),
        ("tracks", "09_tracks", enabled_modules.get("tracks", True)),
    ]
    module_rows = [{"module": label, "status": _receipt_status(results / directory) if enabled else "DISABLED"}
                   for label, directory, enabled in modules]
    summary = outdir / "run_summary.tsv"
    with summary.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["module", "status"], delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(module_rows)

    contrast_rows = _contrast_summary(plan, results, contrasts)
    contrast_summary = outdir / "contrast_summary.tsv"
    _write_tsv(contrast_summary, contrast_rows, CONTRAST_SUMMARY_FIELDS)

    markdown = outdir / "report.md"
    lines = [
        f"# rna_ends2tracks report: {plan.project['project_id']}", "",
        f"- Genomes: `{', '.join(sorted(plan.references))}`",
        f"- Biological samples: {len(plan.samples)}", f"- Pairwise contrasts: {len(contrasts)}",
        "- UMI processing: disabled", "- Coordinate deduplication: disabled", "",
        "## Module status", "", "| Module | Status |", "|---|---|",
    ]
    lines.extend(f"| {row['module']} | {row['status']} |" for row in module_rows)
    lines += ["", "## Contrasts", "",
              "| Contrast | Genome | Numerator | Denominator | Mode | Design |",
              "|---|---|---|---|---|---|"]
    lines.extend(
        f"| {row['contrast_id']} | {row.get('genome', '')} | {row['numerator']} | {row['denominator']} | "
        f"{row.get('design_mode', '')} | `{row.get('resolved_design', '')}` |" for row in contrasts
    )
    lines += [
        "", "## Differential and APA summary", "",
        "| Contrast | DGE significant | DGE up | DGE down | APA-A significant sites | Distal genes | Proximal genes | PCPA |",
        "|---|---:|---:|---:|---:|---:|---:|---:|",
    ]
    lines.extend(
        f"| {row['contrast_id']} | {row['dge_significant']} | {row['dge_up']} | {row['dge_down']} | "
        f"{row['apa_a_significant_sites']} | {row['apa_a_distal_genes']} | "
        f"{row['apa_a_proximal_genes']} | {row['apa_a_pcpa']} |"
        for row in contrast_rows
    )
    lines += [
        "", "## Scientific interpretation", "",
        (
            "Intragenic PAS are candidate premature cleavage/polyadenylation events consistent with premature "
            "transcription termination. QuantSeq alone does not prove loss of downstream engaged RNA polymerase."
        ),
        "", "## Warnings", "",
    ]
    lines.extend(f"- {row.get('warning_code', 'WARNING')}: {row.get('message', '')}" for row in warnings)
    if not warnings:
        lines.append("- No configuration-level warnings.")
    markdown.write_text("\n".join(lines) + "\n", encoding="utf-8")

    c_legend = [
        {"stage": "C0", "meaning": "mapped primary NH=1 alignments", "use": "final BAMs and all-read tracks"},
        {"stage": "C1", "meaning": "exact transcript-end counts", "use": "raw exact-end signal"},
        {"stage": "C1S", "meaning": "uncertain clipped-end counts", "use": "QC only; excluded from PAS calling"},
        {"stage": "C2", "meaning": "filtered exact-end counts", "use": "active-PAS discovery"},
        {"stage": "C2R", "meaning": "internal-priming rejects", "use": "diagnostic tracks and QC"},
        {"stage": "C3", "meaning": "active-PAS counts", "use": "APA and PAS-usage matrices"},
        {"stage": "C4", "meaning": "uniquely assigned active-PAS gene sums", "use": "primary DGE and normalization"},
        {"stage": "C5", "meaning": "conventional exon-overlap counts", "use": "diagnostic comparison only"},
    ]
    warning_rows = warnings or [{"warning_code": "NONE", "message": "No configuration-level warnings."}]
    overview = [
        {"property": "Project", "value": plan.project["project_id"]},
        {"property": "Genomes", "value": ", ".join(sorted(plan.references))},
        {"property": "Biological samples", "value": len(plan.samples)},
        {"property": "Pairwise contrasts", "value": len(contrasts)},
        {"property": "UMI processing", "value": "disabled"},
        {"property": "Coordinate deduplication", "value": "disabled"},
    ]
    body = "".join([
        f"<h1>rna_ends2tracks report: {html.escape(plan.project['project_id'])}</h1>",
        "<p class='note'>QuantSeq REV gene-expression and alternative-polyadenylation summary. "
        "Intragenic PAS are candidate premature cleavage/polyadenylation events; this assay alone "
        "does not prove loss of downstream engaged RNA polymerase.</p>",
        "<h2>Run overview</h2>", _html_table(overview, ["property", "value"]),
        "<h2>Count-universe legend</h2>", _html_table(c_legend, ["stage", "meaning", "use"]),
        "<h2>Module status</h2>", _html_table(module_rows, ["module", "status"]),
        "<h2>Contrast-level results</h2>",
        "<p>Blank cells mean that the corresponding optional analysis was disabled or unavailable.</p>",
        _html_table(contrast_rows, CONTRAST_SUMMARY_FIELDS),
        "<h2>Warnings</h2>", _html_table(warning_rows, ["warning_code", "message"]),
        "<h2>Browser assets</h2><p><a href='IGV_session.xml'>IGV session</a> | "
        "<a href='UCSC_trackDb.txt'>UCSC trackDb</a> | "
        "<a href='contrast_summary.tsv'>full contrast summary TSV</a></p>",
    ])
    html_target.write_text(
        "<!doctype html><html><head><meta charset='utf-8'><title>rna_ends2tracks report</title>"
        "<style>body{font-family:system-ui,sans-serif;max-width:1500px;margin:2rem auto;padding:0 1rem;line-height:1.45}"
        "h1,h2{color:#17324d}.note{background:#eef6fb;border-left:4px solid #2673a6;padding:.8rem 1rem}"
        ".table-wrap{overflow-x:auto;margin-bottom:1.5rem}table{border-collapse:collapse;width:100%;font-size:.9rem}"
        "th,td{border:1px solid #ccd5dd;padding:.4rem .55rem;text-align:left;white-space:nowrap}"
        "th{background:#e9eff4;position:sticky;top:0}tbody tr:nth-child(even){background:#f8fafb}"
        "a{color:#145d91}</style></head>"
        f"<body>{body}</body></html>\n", encoding="utf-8")
    browser = _browser_assets(results, outdir)
    outputs = [markdown, html_target, summary, contrast_summary, *browser]
    write_receipt("report", outdir, signature, outputs, ["rna-ends2tracks", "report"])
    event(log_dir, "report", "completed", str(html_target))
    return html_target
