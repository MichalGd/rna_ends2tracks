from __future__ import annotations

import csv
import json
from pathlib import Path

from .config import signature_for
from .external import event
from .receipts import receipt_valid, write_receipt


def _rows(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        return []
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def make_report(results: Path, force: bool = False) -> Path:
    outdir = results / "08_reports"
    outdir.mkdir(parents=True, exist_ok=True)
    inputs = [results / "00_metadata" / "contrasts.tsv", results / "00_metadata" / "validated_samples.tsv",
              results / "00_metadata" / "resolved_config.json"]
    optional_summary = results / "06_apa_comparison" / "summary.tsv"
    if optional_summary.is_file():
        inputs.append(optional_summary)
    signature = signature_for(inputs, {"module": "report"})
    target = outdir / "report.md"
    if not force and receipt_valid(outdir, signature):
        event(results / "provenance" / "logs", "report", "skipped", "Valid matching receipt")
        return target
    contrasts = _rows(results / "00_metadata" / "contrasts.tsv")
    samples = _rows(results / "00_metadata" / "validated_samples.tsv")
    apa_summary = _rows(results / "06_apa_comparison" / "summary.tsv")
    config = json.loads((results / "00_metadata" / "resolved_config.json").read_text(encoding="utf-8"))
    warnings = [f"{row['contrast_id']}: LOW_REPLICATION_N2 (exploratory)" for row in contrasts if row["design_status"] == "LOW_REPLICATION_N2"]
    warnings.extend(row["message"] for row in _rows(results / "00_metadata" / "warnings.tsv"))
    lines = [
        f"# rna_ends2tracks report: {config['project_id']}", "",
        f"- Species/assembly: `{config['reference']['species']}` / `{config['reference']['assembly']}`",
        f"- Protocol: `{config['protocol']['profile']}`", f"- Biological samples: {len(samples)}",
        f"- Pairwise contrasts: {len(contrasts)}", "- UMI processing: disabled", "- Coordinate deduplication: disabled", "",
        "## Contrasts", "", "| Contrast | Numerator | Denominator | Replicates | Status |", "|---|---:|---:|---:|---|",
    ]
    for row in contrasts:
        lines.append(f"| {row['contrast_id']} | {row['numerator']} | {row['denominator']} | {row['n_num']} / {row['n_den']} | {row['design_status']} |")
    lines += ["", "## APA method reconciliation", ""]
    if apa_summary:
        for row in apa_summary:
            lines.append(f"- {row['class']}: {row['site_count']}")
    else:
        lines.append("APA comparison has not been completed.")
    lines += ["", "## Scientific interpretation", "",
              "Intragenic calls are candidate premature cleavage and polyadenylation events consistent with premature transcription termination. QuantSeq alone does not prove downstream loss of engaged RNA polymerase.", "",
              "## Warnings", ""]
    lines.extend(f"- {warning}" for warning in warnings)
    if not warnings:
        lines.append("- No metadata-level warnings.")
    lines += ["", "## Output index", "",
              "- `03_gene_expression/deseq2/result_index.tsv`", "- `04_apa_a_repository/dexseq/result_index.tsv`",
              "- `05_apa_b_polyaseqtrap_drimseq/drimseq/result_index.tsv`", "- `06_apa_comparison/site_crosswalk.tsv`", ""]
    target.write_text("\n".join(lines), encoding="utf-8")
    write_receipt("report", outdir, signature, [target], ["rna-ends2tracks", "report"])
    event(results / "provenance" / "logs", "report", "completed", str(target))
    return target
