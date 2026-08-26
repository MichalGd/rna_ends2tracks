from __future__ import annotations

import csv
import html
import json
from pathlib import Path

from .config import RunPlan, signature_for, workflow_requirements
from .external import event
from .receipts import receipt_valid, write_receipt


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

    body = "\n".join(f"<p>{html.escape(line)}</p>" for line in lines if line and not line.startswith("|"))
    html_target.write_text(
        "<!doctype html><html><head><meta charset='utf-8'><title>rna_ends2tracks report</title>"
        "<style>body{font-family:sans-serif;max-width:1100px;margin:2rem auto;line-height:1.45}</style></head>"
        f"<body>{body}</body></html>\n", encoding="utf-8")
    browser = _browser_assets(results, outdir)
    outputs = [markdown, html_target, summary, *browser]
    write_receipt("report", outdir, signature, outputs, ["rna-ends2tracks", "report"])
    event(log_dir, "report", "completed", str(html_target))
    return html_target
