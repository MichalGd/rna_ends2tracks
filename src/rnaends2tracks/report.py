from __future__ import annotations

import csv
import html
import json
import re
import shlex
from collections import Counter
from pathlib import Path
from typing import Any

from .config import RunPlan, signature_for, workflow_requirements
from .external import event
from .provenance import generate_provenance_dashboard
from .receipts import receipt_valid, write_receipt


CONTRAST_SUMMARY_FIELDS = [
    "genome", "contrast_id", "numerator", "denominator", "design_mode",
    "paired", "n_pairs", "resolved_design",
    "dge_tested_genes", "dge_significant", "dge_up", "dge_down",
    "apa_a_tested_genes", "apa_a_significant_genes", "apa_a_tested_sites",
    "apa_a_significant_sites", "apa_a_distal_genes",
    "apa_a_proximal_genes", "apa_a_pcpa", "apa_b_tested_sites",
    "apa_b_tested_genes", "apa_b_significant_genes", "apa_b_confirmed_sites",
    "apa_b_distal_genes", "apa_b_proximal_genes", "apa_b_pcpa", "apa_compared_sites",
    "apa_direction_agree", "apa_direction_disagree", "apa_direction_agreement_pct",
]

DGE_EVENT_FIELDS = [
    "genome", "contrast_id", "gene_id", "baseMean", "log2FoldChange",
    "pvalue", "padj", "direction",
]
APA_EVENT_FIELDS = [
    "method", "genome", "contrast_id", "gene_id", "gene_padj", "shift",
    "max_abs_delta_PAU", "weighted_transcript_position_shift_nt",
    "significant_or_confirmed_sites", "pcpa_candidate",
]
ENRICHMENT_TERM_FIELDS = [
    "analysis_type", "genome", "contrast_id", "method", "query", "database",
    "term_id", "term_name", "effect", "padj",
]

TRACK_COLLECTION_ORDER = [
    "all_reads/raw", "all_reads/cpm", "exact_ends/raw", "exact_ends/cpm",
    "filtered_ends/raw", "filtered_ends/cpm", "filtered_ends/deseq2",
    "filtered_ends/robust_cpm", "rejected_ends/raw", "rejected_ends/cpm",
    "active_pas/raw", "active_pas/cpm", "active_pas/deseq2",
    "active_pas/robust_cpm",
]

TRACK_COLLECTION_DESCRIPTIONS = {
    "all_reads/raw": "C0 uniquely mapped read coverage, raw depth",
    "all_reads/cpm": "C0 uniquely mapped read coverage, CPM normalized",
    "exact_ends/raw": "C1 exact transcript 3-prime ends, raw counts",
    "exact_ends/cpm": "C1 exact transcript 3-prime ends, CPM normalized",
    "filtered_ends/raw": "C2 filtered exact ends, raw counts",
    "filtered_ends/cpm": "C2 filtered exact ends, CPM normalized",
    "filtered_ends/deseq2": "C2 filtered exact ends, DESeq2 size-factor normalized",
    "filtered_ends/robust_cpm": "C2 filtered exact ends, DESeq2 robust CPM",
    "rejected_ends/raw": "C2R internal-priming rejects, raw counts",
    "rejected_ends/cpm": "C2R internal-priming rejects, CPM normalized",
    "active_pas/raw": "C3 active polyadenylation sites, raw counts",
    "active_pas/cpm": "C3 active polyadenylation sites, CPM normalized",
    "active_pas/deseq2": "C3 active polyadenylation sites, DESeq2 size-factor normalized",
    "active_pas/robust_cpm": "C3 active polyadenylation sites, DESeq2 robust CPM",
}

TRACK_COLLECTION_COLORS = [
    "0,102,204", "0,153,153", "76,114,176", "85,168,104",
    "196,78,82", "237,201,72", "175,122,161", "255,157,167",
    "156,117,95", "186,176,172", "31,119,180", "44,160,44",
    "214,39,40", "148,103,189",
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
                "paired": contrast.get("paired", ""),
                "n_pairs": contrast.get("n_pairs", ""),
                "resolved_design": contrast.get("resolved_design", ""),
                "dge_tested_genes": "",
                "dge_significant": "",
                "dge_up": "",
                "dge_down": "",
                "apa_a_tested_sites": "",
                "apa_a_tested_genes": "",
                "apa_a_significant_genes": "",
                "apa_a_significant_sites": "",
                "apa_a_distal_genes": "",
                "apa_a_proximal_genes": "",
                "apa_a_pcpa": "" if apa_a_pcpa is None else apa_a_pcpa[contrast_id],
                "apa_b_tested_sites": "",
                "apa_b_tested_genes": "",
                "apa_b_significant_genes": "",
                "apa_b_confirmed_sites": "",
                "apa_b_distal_genes": "",
                "apa_b_proximal_genes": "",
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
                gene_summary = Path(index.get("gene_summary_file", ""))
                if gene_summary.is_file():
                    genes = _required_rows(gene_summary, "APA-A gene summary")
                    row.update({
                        "apa_a_tested_genes": len(genes),
                        "apa_a_significant_genes": sum(
                            _number(item.get("gene_padj")) is not None
                            and _number(item.get("gene_padj")) <= fdr
                            for item in genes
                        ),
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
                gene_summary = Path(index.get("gene_summary_file", ""))
                if gene_summary.is_file():
                    genes = _required_rows(gene_summary, "APA-B gene summary")
                    significant_genes = [
                        item for item in genes
                        if _number(item.get("gene_padj")) is not None
                        and _number(item.get("gene_padj")) <= fdr
                    ]
                    shifts = Counter(item.get("shift", "") for item in significant_genes)
                    row.update({
                        "apa_b_tested_genes": len(genes),
                        "apa_b_significant_genes": len(significant_genes),
                        "apa_b_distal_genes": shifts["distal"],
                        "apa_b_proximal_genes": shifts["proximal"],
                    })
            row.update(comparison.get(contrast_id, {}))
            compared = int(row.get("apa_compared_sites") or 0)
            agreed = int(row.get("apa_direction_agree") or 0)
            row["apa_direction_agreement_pct"] = (
                f"{100 * agreed / compared:.2f}" if compared else ""
            )
            summary.append(row)
    return summary


def _write_tsv(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(rows)


def _html_table(
    rows: list[dict[str, Any]], fields: list[str], table_id: str = "",
) -> str:
    headings = "".join(f"<th>{html.escape(field.replace('_', ' '))}</th>" for field in fields)
    body = "".join(
        "<tr>" + "".join(f"<td>{html.escape(str(row.get(field, '')))}</td>" for field in fields) + "</tr>"
        for row in rows
    )
    identifier = f" id='{html.escape(table_id)}'" if table_id else ""
    return f"<div class='table-wrap'><table{identifier}><thead><tr>{headings}</tr></thead><tbody>{body}</tbody></table></div>"


def _image_gallery(outdir: Path, images: list[tuple[str, Path]]) -> str:
    cards: list[str] = []
    for label, path in images:
        if not path.is_file():
            continue
        href = path.relative_to(outdir).as_posix() if path.is_relative_to(outdir) else "../" + path.relative_to(outdir.parent).as_posix()
        cards.append(
            "<figure><a href='" + html.escape(href) + "'><img loading='lazy' src='" + html.escape(href) +
            "' alt='" + html.escape(label) + "'></a><figcaption>" + html.escape(label) + "</figcaption></figure>"
        )
    return "<div class='gallery'>" + "".join(cards) + "</div>" if cards else "<p>No plot was available.</p>"


def _apa_b_interpretation(plan: RunPlan) -> tuple[str, list[dict[str, Any]]]:
    settings = plan.project.get("apa_b", {})
    if not settings.get("enabled", False):
        return "DISABLED_NOT_VALIDATED", [{
            "property": "Interpretation status", "value": "DISABLED_NOT_VALIDATED",
        }, {
            "property": "Meaning",
            "value": "APA-B results were not generated and must not be inferred from APA-A.",
        }]
    manifest = Path(str(settings.get("validation_manifest", "")))
    try:
        value = json.loads(manifest.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return "INVALID_VALIDATION_MANIFEST", [{"property": "Validation manifest", "value": str(manifest)}]
    engine = value.get("engine", {})
    return "VALIDATED_PILOT_ACCEPTED", [
        {"property": "Interpretation status", "value": "VALIDATED_PILOT_ACCEPTED"},
        {"property": "Engine", "value": engine.get("name", "")},
        {"property": "Pinned source commit", "value": engine.get("source_commit", "")},
        {"property": "Reviewed by", "value": value.get("reviewed_by", "")},
        {"property": "Accepted at", "value": value.get("accepted_at", "")},
        {"property": "Assemblies", "value": ", ".join(value.get("assemblies", []))},
        {"property": "UMI / coordinate deduplication", "value": "disabled / disabled"},
    ]


def _apa_b_gene_events(results: Path, fdr: float, limit: int = 100) -> list[dict[str, Any]]:
    events: list[dict[str, Any]] = []
    root = results / "07_apa_b"
    for index_path in sorted(root.glob("*/drimseq/result_index.tsv")) if root.is_dir() else []:
        genome = index_path.relative_to(root).parts[0]
        for index in _rows(index_path):
            summary = Path(index.get("gene_summary_file", ""))
            if not summary.is_file():
                continue
            for row in _rows(summary):
                adjusted = _number(row.get("gene_padj"))
                if adjusted is None or adjusted > fdr:
                    continue
                events.append({
                    "genome": genome, "contrast_id": index.get("contrast_id", ""),
                    "gene_id": row.get("gene_id", ""), "gene_padj": adjusted,
                    "shift": row.get("shift", ""),
                    "max_abs_delta_PAU": row.get("max_abs_delta_PAU", ""),
                    "weighted_transcript_position_shift_nt": row.get("weighted_transcript_position_shift_nt", ""),
                    "confirmed_sites": row.get("confirmed_sites", ""),
                })
    events.sort(key=lambda row: float(row["gene_padj"]))
    return events[:limit]


def _top_dge_events(results: Path, fdr: float, limit_per_contrast: int = 25) -> list[dict[str, Any]]:
    events: list[dict[str, Any]] = []
    root = results / "05_gene_expression"
    # C4 active-PAS gene sums are the primary DGE result.  C5 featureCounts is
    # diagnostic only and must not be mixed into the biological summary.
    for index_path in sorted(root.glob("*/C4_primary_deseq2/result_index.tsv")) if root.is_dir() else []:
        genome = index_path.relative_to(root).parts[0]
        for index in _rows(index_path):
            result = Path(index.get("result_file", ""))
            if not result.is_file():
                continue
            selected = []
            for row in _rows(result):
                padj = _number(row.get("padj"))
                effect = _number(row.get("log2FoldChange"))
                if padj is None or padj > fdr or effect is None:
                    continue
                selected.append((padj, -abs(effect), row.get("gene_id", ""), row, effect))
            for _padj, _rank, _gene, row, effect in sorted(selected, key=lambda item: item[:3])[:limit_per_contrast]:
                events.append({
                    "genome": genome, "contrast_id": index.get("contrast_id", ""),
                    "gene_id": row.get("gene_id", ""), "baseMean": row.get("baseMean", ""),
                    "log2FoldChange": row.get("log2FoldChange", ""),
                    "pvalue": row.get("pvalue", ""), "padj": row.get("padj", ""),
                    "direction": "up" if effect > 0 else "down",
                })
    return events


def _top_apa_events(results: Path, fdr: float, limit_per_contrast: int = 25) -> list[dict[str, Any]]:
    events: list[dict[str, Any]] = []
    definitions = (
        ("APA-A", results / "06_apa_a_mcell2019", "significant_sites"),
        ("APA-B", results / "07_apa_b", "confirmed_sites"),
    )
    for method, root, site_field in definitions:
        if not root.is_dir():
            continue
        pcpa_by_genome: dict[str, set[tuple[str, str]]] = {}
        for pcpa_path in root.glob("*/candidate_pcpa.tsv"):
            genome = pcpa_path.parent.name
            pcpa_by_genome[genome] = {
                (row.get("contrast_id", ""), row.get("gene_id", "")) for row in _rows(pcpa_path)
            }
        for index_path in sorted(root.rglob("result_index.tsv")):
            genome = index_path.relative_to(root).parts[0]
            for index in _rows(index_path):
                summary = Path(index.get("gene_summary_file", ""))
                if not summary.is_file():
                    continue
                selected = []
                for row in _rows(summary):
                    padj = _number(row.get("gene_padj"))
                    effect = _number(row.get("max_abs_delta_PAU")) or 0.0
                    if padj is not None and padj <= fdr:
                        selected.append((padj, -abs(effect), row.get("gene_id", ""), row))
                for _padj, _effect, _gene, row in sorted(
                    selected, key=lambda item: item[:3]
                )[:limit_per_contrast]:
                    contrast_id = index.get("contrast_id", "")
                    gene_id = row.get("gene_id", "")
                    events.append({
                        "method": method, "genome": genome, "contrast_id": contrast_id,
                        "gene_id": gene_id, "gene_padj": row.get("gene_padj", ""),
                        "shift": row.get("shift", ""),
                        "max_abs_delta_PAU": row.get("max_abs_delta_PAU", ""),
                        "weighted_transcript_position_shift_nt": row.get(
                            "weighted_transcript_position_shift_nt", ""),
                        "significant_or_confirmed_sites": row.get(site_field, ""),
                        "pcpa_candidate": str(
                            (contrast_id, gene_id) in pcpa_by_genome.get(genome, set())
                        ).lower(),
                    })
    return events


def _top_enrichment_terms(
    enrichment_rows: list[dict[str, str]], fdr: float, limit_per_job: int = 10,
) -> list[dict[str, Any]]:
    terms: list[dict[str, Any]] = []
    for index in enrichment_rows:
        for method, field in (("ORA", "ora_file"), ("GSEA", "gsea_file")):
            selected: list[tuple[float, dict[str, Any]]] = []
            source = Path(index.get(field, ""))
            if not source.is_file():
                continue
            for row in _rows(source):
                padj = _number(row.get("padj"))
                if padj is None or padj > fdr:
                    continue
                effect = row.get("NES", "") if method == "GSEA" else row.get("overlap_count", "")
                selected.append((padj, {
                    "analysis_type": index.get("analysis_type", ""),
                    "genome": index.get("genome", ""),
                    "contrast_id": index.get("contrast_id", ""), "method": method,
                    "query": row.get("query", ""), "database": row.get("database", ""),
                    "term_id": row.get("term_id", ""), "term_name": row.get("term_name", ""),
                    "effect": effect, "padj": row.get("padj", ""),
                }))
            # Keep both analyses visible: a large ORA result must not crowd all
            # ranked GSEA terms (or vice versa) out of the final report.
            terms.extend(
                row for _padj, row in sorted(selected, key=lambda item: item[0])[:limit_per_job]
            )
    return terms


def _track_collections(results: Path) -> tuple[list[Path], list[dict[str, Any]]]:
    track_root = results / "09_tracks"
    bigwigs = sorted(track_root.rglob("*.bw")) if track_root.is_dir() else []
    grouped: dict[str, list[Path]] = {}
    for path in bigwigs:
        relative = path.relative_to(track_root)
        group = "/".join(relative.parts[:-1]) or "ungrouped"
        grouped.setdefault(group, []).append(path)
    order = [group for group in TRACK_COLLECTION_ORDER if group in grouped]
    order.extend(sorted(set(grouped).difference(order)))
    rows = [{
        "collection": group,
        "description": TRACK_COLLECTION_DESCRIPTIONS.get(group, "Generated BigWig collection"),
        "bigwigs": len(grouped[group]),
        "transcript_plus": sum("transcript_plus" in path.name for path in grouped[group]),
        "transcript_minus": sum("transcript_minus" in path.name for path in grouped[group]),
    } for group in order]
    return bigwigs, rows


def _ucsc_track_line(
    path: Path, track_root: Path, index: int, color: str, description: str,
    public_prefix: str, relative_prefix: str, negate_minus: bool, view_limits: str,
) -> str:
    relative = path.relative_to(track_root).as_posix()
    # The optional public destination is deliberately flat: generated BigWig
    # basenames already contain sample, family, normalization and strand.
    url = f"{public_prefix}/{path.name}" if public_prefix else f"{relative_prefix}/{relative}"
    # UCSC custom-track labels are intentionally conservative: name <=15 and
    # description <=60 characters. The full filename remains in bigDataUrl.
    clean_description = f"{description}: {path.name}".replace('"', "'")[:60]
    fields = [
        "track", "type=bigWig", f'name="rna_ends_{index:04d}"',
        f'description="{clean_description}"', f"bigDataUrl={url}", "visibility=full",
        f"color={color}", f"priority={index}", "autoScale=off", "alwaysZero=on",
        "gridDefault=on", "graphType=bar", "windowingFunction=mean",
        f"viewLimits={view_limits}", "maxHeightPixels=100:50:8",
    ]
    if negate_minus and "transcript_minus" in path.name:
        fields.append("negateValues=on")
    return " ".join(fields)


def _validate_ucsc_track_lines(lines: list[str], source: str = "descriptor") -> int:
    """Validate the one-line UCSC bigWig custom-track contract."""
    required = {"type", "name", "description", "bigDataUrl"}
    names: set[str] = set()
    count = 0
    for number, raw in enumerate(lines, 1):
        if "\n" in raw or "\r" in raw:
            raise RuntimeError(f"UCSC descriptor is not one line at {source}:{number}")
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        try:
            tokens = shlex.split(line, posix=True)
        except ValueError as exc:
            raise RuntimeError(f"Invalid UCSC quoting at {source}:{number}: {exc}") from exc
        if not tokens or tokens[0] != "track":
            raise RuntimeError(f"UCSC descriptor must start with 'track' at {source}:{number}")
        attributes: dict[str, str] = {}
        for token in tokens[1:]:
            if "=" not in token:
                raise RuntimeError(f"Invalid UCSC key=value token at {source}:{number}: {token}")
            key, value = token.split("=", 1)
            if not key or not value or key in attributes:
                raise RuntimeError(f"Invalid or duplicate UCSC attribute at {source}:{number}: {key}")
            attributes[key] = value
        missing = sorted(required.difference(attributes))
        if missing:
            raise RuntimeError(f"Missing UCSC attributes at {source}:{number}: {', '.join(missing)}")
        if attributes["type"] != "bigWig":
            raise RuntimeError(f"UCSC track type must be bigWig at {source}:{number}")
        name = attributes["name"]
        if len(name) > 15 or name in names:
            raise RuntimeError(f"Invalid or duplicate UCSC track name at {source}:{number}: {name}")
        names.add(name)
        if len(attributes["description"]) > 60:
            raise RuntimeError(f"UCSC description exceeds 60 characters at {source}:{number}")
        url = attributes["bigDataUrl"]
        if (any(character.isspace() for character in url)
                or any(character in url for character in "[]")
                or not (url.startswith(("http://", "https://", "../", "./")))
                or not url.lower().endswith((".bw", ".bigwig"))):
            raise RuntimeError(f"Invalid UCSC bigDataUrl at {source}:{number}: {url}")
        color = attributes.get("color", "")
        if color:
            components = color.split(",")
            if len(components) != 3 or any(
                not component.isdigit() or not 0 <= int(component) <= 255
                for component in components
            ):
                raise RuntimeError(f"Invalid UCSC RGB color at {source}:{number}: {color}")
        view_limits = attributes.get("viewLimits", "")
        if view_limits and not re.fullmatch(
            r"-?[0-9]+(?:\.[0-9]+)?:-?[0-9]+(?:\.[0-9]+)?", view_limits,
        ):
            raise RuntimeError(f"Invalid UCSC viewLimits at {source}:{number}: {view_limits}")
        if attributes.get("negateValues") not in {None, "on", "off"}:
            raise RuntimeError(f"Invalid UCSC negateValues at {source}:{number}")
        for key in ("autoScale", "alwaysZero", "gridDefault"):
            if attributes.get(key) not in {None, "on", "off"}:
                raise RuntimeError(f"Invalid UCSC {key} at {source}:{number}")
        if attributes.get("priority") and not re.fullmatch(r"[0-9]+(?:\.[0-9]+)?", attributes["priority"]):
            raise RuntimeError(f"Invalid UCSC priority at {source}:{number}")
        if attributes.get("maxHeightPixels") and not re.fullmatch(
            r"[0-9]+:[0-9]+:[0-9]+", attributes["maxHeightPixels"],
        ):
            raise RuntimeError(f"Invalid UCSC maxHeightPixels at {source}:{number}")
        count += 1
    if not count:
        raise RuntimeError(f"No UCSC track definitions were generated in {source}")
    return count


def _browser_assets(
    plan: RunPlan, results: Path, outdir: Path,
) -> tuple[list[Path], list[dict[str, Any]]]:
    bigwigs, collection_rows = _track_collections(results)
    track_root = results / "09_tracks"
    settings = plan.project.get("tracks", {})
    public_prefix = str(settings.get("ucsc_bigdata_url_prefix", "")).rstrip("/")
    negate_minus = bool(settings.get("ucsc_negate_minus_tracks", True))
    view_limits = str(settings.get("ucsc_view_limits", "0:12"))

    inventory = outdir / "bigwig_collections.txt"
    with inventory.open("w", encoding="utf-8") as handle:
        for row in collection_rows:
            group = str(row["collection"])
            handle.write(f"[{group}]\n")
            for path in bigwigs:
                if path.relative_to(track_root).parent.as_posix() == group:
                    handle.write(path.name + "\n")
            handle.write("\n")

    descriptor_dir = outdir / "ucsc_track_descriptors"
    descriptor_dir.mkdir(parents=True, exist_ok=True)
    # Alpha.10 originally placed this compatibility file beside the folder.
    # Remove that workflow-owned stale copy so every UCSC descriptor has one home.
    (outdir / "UCSC_trackDb.txt").unlink(missing_ok=True)
    # Remove only filenames owned by previous exporter layouts. This prevents
    # stale collection/family descriptors surviving a forced report rebuild.
    known_groups = set(TRACK_COLLECTION_ORDER) | {
        str(row["collection"]) for row in collection_rows
    }
    known_families = {group.split("/", 1)[0] for group in known_groups}
    for stale in [
        *(descriptor_dir / f"{group.replace('/', '__')}.txt" for group in known_groups),
        *(descriptor_dir / f"{family}.txt" for family in known_families),
    ]:
        stale.unlink(missing_ok=True)
    combined = descriptor_dir / "UCSC_bigWig_tracks.oneline.txt"
    generated: list[Path] = [inventory, combined]
    all_lines: list[str] = []
    family_lines: dict[str, list[str]] = {}
    track_index = 0
    for group_index, row in enumerate(collection_rows):
        group = str(row["collection"])
        family = group.split("/", 1)[0]
        description = str(row["description"])
        color = TRACK_COLLECTION_COLORS[group_index % len(TRACK_COLLECTION_COLORS)]
        group_lines = [f"# {group} - {description}"]
        for path in bigwigs:
            if path.relative_to(track_root).parent.as_posix() != group:
                continue
            track_index += 1
            group_lines.append(_ucsc_track_line(
                path, track_root, track_index, color, description, public_prefix,
                "../../09_tracks", negate_minus, view_limits,
            ))
        family_lines.setdefault(family, []).extend([*group_lines, ""])
        all_lines.extend([*group_lines, ""])
    if all_lines:
        _validate_ucsc_track_lines(all_lines, combined.name)
        combined.write_text("\n".join(all_lines).rstrip() + "\n", encoding="utf-8")
    else:
        combined.write_text("# No BigWig tracks were generated.\n", encoding="utf-8")

    validation_rows: list[dict[str, Any]] = []
    for family, lines in family_lines.items():
        target = descriptor_dir / f"{family}.txt"
        count = _validate_ucsc_track_lines(lines, target.name)
        target.write_text("\n".join(lines).rstrip() + "\n", encoding="utf-8")
        generated.append(target)
        validation_rows.append({"file": target.name, "track_lines": count, "status": "PASS"})
    trackdb = descriptor_dir / "UCSC_trackDb.txt"
    trackdb.write_text(combined.read_text(encoding="utf-8"), encoding="utf-8")
    generated.append(trackdb)
    combined_status = "PASS" if track_index else "NO_TRACKS"
    validation_rows.extend([
        {"file": combined.name, "track_lines": track_index, "status": combined_status},
        {"file": trackdb.name, "track_lines": track_index, "status": combined_status},
    ])
    validation = descriptor_dir / "UCSC_descriptor_validation.tsv"
    _write_tsv(validation, validation_rows, ["file", "track_lines", "status"])
    generated.append(validation)
    igv = outdir / "IGV_session.xml"
    resources = "\n".join(
        f'    <Resource path="../{html.escape(path.relative_to(results).as_posix())}"/>' for path in bigwigs
    )
    igv.write_text(f'<?xml version="1.0" encoding="UTF-8"?>\n<Session version="8">\n  <Resources>\n{resources}\n  </Resources>\n</Session>\n',
                   encoding="utf-8")
    generated.append(igv)
    return generated, collection_rows


def _main_artifacts(results: Path, outdir: Path) -> list[dict[str, str]]:
    candidates: list[tuple[str, Path]] = [
        ("MultiQC report", results / "01_qc" / "multiqc" / "multiqc_report.html"),
        ("RSeQC MultiQC report", results / "01_qc" / "rseqc" / "multiqc" / "multiqc_report.html"),
        ("RSeQC sample summary", results / "01_qc" / "rseqc" / "rseqc_summary.tsv"),
        ("RSeQC gene-body coverage", results / "01_qc" / "rseqc" / "gene_body_coverage.svg"),
        ("Protocol-orientation audit", results / "02_alignment" / "protocol_orientation.tsv"),
        ("Track normalization table", results / "09_tracks" / "track_normalization.tsv"),
        ("One-column BigWig inventory", outdir / "bigwig_collections.txt"),
        ("Combined one-line UCSC descriptors", outdir / "ucsc_track_descriptors" / "UCSC_bigWig_tracks.oneline.txt"),
        ("UCSC descriptor validation", outdir / "ucsc_track_descriptors" / "UCSC_descriptor_validation.tsv"),
        ("IGV session", outdir / "IGV_session.xml"),
        ("Full contrast summary", outdir / "contrast_summary.tsv"),
        ("Differential gene-expression summary", outdir / "differential_gene_expression_summary.tsv"),
        ("Top differential genes", outdir / "top_differential_genes.tsv"),
        ("Alternative-polyadenylation summary", outdir / "alternative_polyadenylation_summary.tsv"),
        ("Top APA gene events", outdir / "top_apa_gene_events.tsv"),
        ("Top enrichment terms", outdir / "top_enrichment_terms.tsv"),
        ("Enrichment result index", outdir / "enrichment_summary" / "enrichment_index.tsv"),
        ("Provenance dashboard", outdir / "provenance_dashboard" / "dashboard.json"),
        ("Receipt inventory", outdir / "provenance_dashboard" / "receipt_inventory.tsv"),
        ("Environment packages", outdir / "provenance_dashboard" / "environment_packages.tsv"),
        ("Complete output manifest", outdir / "provenance_dashboard" / "output_manifest.tsv"),
    ]
    for root, label in (
        (results / "05_gene_expression", "DGE result index"),
        (results / "06_apa_a_mcell2019", "APA-A result index"),
        (results / "07_apa_b", "APA-B result index"),
    ):
        for path in (sorted(root.rglob("result_index.tsv")) if root.is_dir() else []):
            genome = path.relative_to(root).parts[0]
            candidates.append((f"{genome} {label}", path))
    rows = []
    for label, path in candidates:
        if not path.is_file():
            continue
        if path.is_relative_to(outdir):
            href = path.relative_to(outdir).as_posix()
        else:
            href = "../" + path.relative_to(results).as_posix()
        rows.append({"artifact": label, "path": href})
    return rows


def _star_qc_rows(results: Path) -> list[dict[str, str]]:
    root = results / "02_alignment"
    rows: list[dict[str, str]] = []
    if not root.is_dir():
        return rows
    wanted = {
        "Number of input reads": "input_reads",
        "Uniquely mapped reads number": "uniquely_mapped_reads",
        "Uniquely mapped reads %": "uniquely_mapped_pct",
        "% of reads mapped to multiple loci": "multimapped_pct",
        "% of reads unmapped: too short": "unmapped_too_short_pct",
    }
    for path in sorted(root.rglob("*Log.final.out")):
        metrics: dict[str, str] = {}
        for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
            if "|" not in line:
                continue
            key, value = (part.strip() for part in line.split("|", 1))
            if key in wanted:
                metrics[wanted[key]] = value
        relative = path.relative_to(root)
        rows.append({
            "sample_id": relative.parts[0] if len(relative.parts) > 1 else "",
            "lane": path.name.removesuffix(".star.Log.final.out"),
            **{field: metrics.get(field, "") for field in wanted.values()},
        })
    return rows


def _exact_funnel_rows(results: Path) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    root = results / "03_exact_ends"
    for path in (sorted(root.glob("*/*/end_audit.json")) if root.is_dir() else []):
        payload = json.loads(path.read_text(encoding="utf-8"))
        c0 = int(payload.get("C0", 0))
        c1 = int(payload.get("C1", 0))
        c2 = int(payload.get("C2", 0))
        rows.append({
            "sample_id": payload.get("sample_id", path.parent.name),
            "genome": payload.get("genome", path.parent.parent.name),
            "C0": c0, "C1": c1, "C1S": int(payload.get("C1S", 0)),
            "C2": c2, "C2R": int(payload.get("C2R", 0)),
            "C1_over_C0_pct": f"{100 * c1 / c0:.2f}" if c0 else "",
            "C2_over_C1_pct": f"{100 * c2 / c1:.2f}" if c1 else "",
            "mask_rescued": int(payload.get("mask_rescued", 0)),
        })
    return rows


def _active_pas_rows(results: Path) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    root = results / "04_active_pas"
    for path in (sorted(root.glob("*/count_universe_audit.json")) if root.is_dir() else []):
        payload = json.loads(path.read_text(encoding="utf-8"))
        c2 = int(payload.get("C2_total", 0))
        c3 = int(payload.get("C3_total", 0))
        rows.append({
            "genome": payload.get("genome", path.parent.name),
            "samples": payload.get("samples", ""),
            "active_pas": payload.get("active_pas", ""),
            "C2_total": c2, "C3_total": c3,
            "C3_over_C2_pct": f"{100 * c3 / c2:.2f}" if c2 else "",
            "ambiguous_pas": payload.get("ambiguous_pas", ""),
            "pcpa_candidate_sites": payload.get("pcpa_candidate_sites", ""),
        })
    return rows


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
    inputs.extend(Path(str(plan.project.get(key, ""))) for key in ("_config_path", "_samplesheet_path"))
    inputs.extend(
        results / directory / "run_receipt.json"
        for directory in (
            "02_alignment", "03_exact_ends", "04_active_pas", "05_gene_expression",
            "06_apa_a_mcell2019", "07_apa_b", "08_apa_comparison", "09_tracks",
            "01_qc/rseqc",
        )
    )
    result_indexes = sorted((results / "05_gene_expression").rglob("result_index.tsv"))
    result_indexes.extend(sorted((results / "06_apa_a_mcell2019").rglob("result_index.tsv")))
    result_indexes.extend(sorted((results / "07_apa_b").rglob("result_index.tsv")))
    inputs.extend(result_indexes)
    for index_path in result_indexes:
        for row in _rows(index_path):
            for field in ("result_file", "shift_file", "gene_summary_file"):
                if row.get(field):
                    inputs.append(Path(row[field]))
    inputs.extend(sorted((results / "06_apa_a_mcell2019").rglob("candidate_pcpa.tsv")))
    inputs.extend(sorted((results / "07_apa_b").rglob("candidate_pcpa.tsv")))
    inputs.append(results / "08_apa_comparison" / "effect_concordance.tsv")
    inputs.extend(sorted((results / "02_alignment").rglob("*Log.final.out")))
    inputs.append(results / "02_alignment" / "protocol_orientation.tsv")
    inputs.extend([
        results / "01_qc" / "fastq_screen" / "fastq_screen_summary.tsv",
        results / "01_qc" / "fastq_screen" / "fastq_screen_metrics.tsv",
        results / "01_qc" / "rseqc" / "rseqc_summary.tsv",
        results / "01_qc" / "rseqc" / "gene_body_coverage.tsv",
        results / "01_qc" / "rseqc" / "gene_body_coverage.svg",
    ])
    inputs.extend(sorted((results / "03_exact_ends").glob("*/*/end_audit.json")))
    inputs.extend(sorted((results / "04_active_pas").glob("*/count_universe_audit.json")))
    inputs.extend(sorted((results / "05_gene_expression").rglob("*.png")))
    enrichment_index_path = results / "10_reports" / "enrichment_summary" / "enrichment_index.tsv"
    inputs.extend(sorted((results / "10_reports" / "enrichment_summary").rglob("*.tsv")))
    for row in _rows(enrichment_index_path):
        for field in (
            "prepared_gene_table", "ora_file", "gsea_file", "plot_index", "mapping_audit", "provenance_file",
        ):
            if row.get(field):
                inputs.append(Path(row[field]))
        plot_index = Path(row.get("plot_index", ""))
        if plot_index.is_file():
            for plot in _rows(plot_index):
                inputs.extend(Path(plot[key]) for key in ("pdf", "png") if plot.get(key))
    inputs.append(results / "10_reports" / "enrichment_summary" / "run_receipt.json")
    validation_manifest = Path(str(plan.project.get("apa_b", {}).get("validation_manifest", "")))
    if validation_manifest.is_file():
        inputs.append(validation_manifest)
    inputs = [path for path in inputs if path.is_file()]
    signature = signature_for(inputs, {
        "module": "report", "project": plan.project["project_id"],
        "reporting": plan.project.get("reporting", {}),
        "browser_tracks": {
            key: plan.project.get("tracks", {}).get(key)
            for key in ("ucsc_bigdata_url_prefix", "ucsc_negate_minus_tracks", "ucsc_view_limits")
        },
        "rseqc": plan.project.get("rseqc", {}),
    })
    html_target = outdir / "report.html"
    if not force and receipt_valid(outdir, signature):
        event(log_dir, "report", "skipped", "Valid matching receipt")
        return html_target

    contrasts = _rows(results / "00_metadata" / "contrasts.tsv")
    warnings = _rows(results / "00_metadata" / "warnings.tsv")
    enabled_modules = plan.project.get("modules", {})
    requirements = workflow_requirements(plan.project)
    enrichment_enabled = bool(
        (enabled_modules.get("dge_enrichment", False) and enabled_modules.get("gene_expression", True))
        or (enabled_modules.get("apa_enrichment", False) and (
            enabled_modules.get("apa_a", True) or plan.project.get("apa_b", {}).get("enabled", False)
        ))
    )
    modules = [
        ("alignment", "02_alignment", True),
        ("RSeQC", "01_qc/rseqc", enabled_modules.get("rseqc", False)),
        ("exact ends", "03_exact_ends", requirements["exact_ends"]),
        ("active PAS", "04_active_pas", requirements["active_pas"]),
        ("gene expression", "05_gene_expression", enabled_modules.get("gene_expression", True)),
        ("APA-A", "06_apa_a_mcell2019", enabled_modules.get("apa_a", True)),
        ("APA-B", "07_apa_b", plan.project.get("apa_b", {}).get("enabled", False)),
        ("APA comparison", "08_apa_comparison", requirements["apa_comparison"]),
        ("gene-set enrichment", "10_reports/enrichment_summary",
         enrichment_enabled),
        ("tracks", "09_tracks", enabled_modules.get("tracks", True)),
    ]
    module_rows = [{"module": label, "status": _receipt_status(results / directory) if enabled else "DISABLED"}
                   for label, directory, enabled in modules]
    fastq_screen_rows = _rows(results / "01_qc" / "fastq_screen" / "fastq_screen_summary.tsv")
    fastq_screen_metrics = _rows(results / "01_qc" / "fastq_screen" / "fastq_screen_metrics.tsv")
    fastq_screen_statuses = {row.get("status", "") for row in fastq_screen_rows}
    fastq_screen_status = (
        "PASS" if fastq_screen_statuses == {"PASS"}
        else "DISABLED" if fastq_screen_statuses == {"DISABLED"}
        else "SKIPPED_MISSING_CONFIG" if "SKIPPED_MISSING_CONFIG" in fastq_screen_statuses
        else "UNAVAILABLE"
    )
    module_rows.insert(1, {"module": "FastQ Screen", "status": fastq_screen_status})
    summary = outdir / "run_summary.tsv"
    with summary.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["module", "status"], delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(module_rows)

    contrast_rows = _contrast_summary(plan, results, contrasts)
    contrast_summary = outdir / "contrast_summary.tsv"
    _write_tsv(contrast_summary, contrast_rows, CONTRAST_SUMMARY_FIELDS)
    dge_summary = outdir / "differential_gene_expression_summary.tsv"
    dge_fields = [
        "genome", "contrast_id", "numerator", "denominator", "design_mode",
        "paired", "n_pairs", "resolved_design",
        "dge_tested_genes", "dge_significant", "dge_up", "dge_down",
    ]
    _write_tsv(dge_summary, contrast_rows, dge_fields)
    apa_summary = outdir / "alternative_polyadenylation_summary.tsv"
    apa_fields = [
        field for field in CONTRAST_SUMMARY_FIELDS
        if field in {
            "genome", "contrast_id", "numerator", "denominator", "design_mode",
            "paired", "n_pairs", "resolved_design",
        }
        or field.startswith("apa_")
    ]
    _write_tsv(apa_summary, contrast_rows, apa_fields)
    browser, track_rows = _browser_assets(plan, results, outdir)
    samples = _rows(results / "00_metadata" / "validated_samples.tsv")
    star_rows = _star_qc_rows(results)
    orientation_rows = _rows(results / "02_alignment" / "protocol_orientation.tsv")
    rseqc_rows = _rows(results / "01_qc" / "rseqc" / "rseqc_summary.tsv")
    funnel_rows = _exact_funnel_rows(results)
    active_pas_rows = _active_pas_rows(results)
    enrichment_rows = _rows(outdir / "enrichment_summary" / "enrichment_index.tsv")
    report_fdr = float(plan.project.get("reporting", {}).get("fdr", 0.05))
    enrichment_fdr = float(plan.project.get("enrichment", {}).get("padj", 0.05))
    dge_events = _top_dge_events(results, report_fdr)
    apa_events = _top_apa_events(results, report_fdr)
    enrichment_terms = _top_enrichment_terms(enrichment_rows, enrichment_fdr)
    dge_events_path = outdir / "top_differential_genes.tsv"
    apa_events_path = outdir / "top_apa_gene_events.tsv"
    enrichment_terms_path = outdir / "top_enrichment_terms.tsv"
    _write_tsv(dge_events_path, dge_events, DGE_EVENT_FIELDS)
    _write_tsv(apa_events_path, apa_events, APA_EVENT_FIELDS)
    _write_tsv(enrichment_terms_path, enrichment_terms, ENRICHMENT_TERM_FIELDS)
    provenance_outputs = generate_provenance_dashboard(plan, results, outdir)
    artifact_rows = _main_artifacts(results, outdir)
    apa_b_status, apa_b_rows = _apa_b_interpretation(plan)
    apa_b_gene_events = _apa_b_gene_events(
        results, float(plan.project.get("reporting", {}).get("fdr", 0.05)),
    ) if apa_b_status == "VALIDATED_PILOT_ACCEPTED" else []
    plot_images: list[tuple[str, Path]] = []
    for genome in sorted(plan.references):
        primary = results / "05_gene_expression" / genome / "C4_primary_deseq2"
        plot_images.extend([
            (f"{genome} C4 variance-stabilized PCA", primary / "C4_vst_pca.png"),
            (f"{genome} sample-distance heatmap", primary / "C4_sample_distances.png"),
        ])
        for row in _rows(primary / "result_index.tsv"):
            for field, label in (("ma_png", "MA"), ("volcano_png", "volcano")):
                if row.get(field):
                    plot_images.append((f"{genome} {row['contrast_id']} {label}", Path(row[field])))
    enrichment_images = [
        (f"{row.get('genome', '')} {row.get('contrast_id', '')} {row.get('analysis_type', '')} enrichment",
         Path(row["plot_png"]))
        for row in enrichment_rows if row.get("plot_png")
    ]
    for row in enrichment_rows:
        plot_index = Path(row.get("plot_index", ""))
        for plot in (_rows(plot_index) if plot_index.is_file() else []):
            if plot.get("png"):
                enrichment_images.append((
                    f"{row.get('analysis_type', '')} {row.get('contrast_id', '')} "
                    f"{plot.get('method', '')} {plot.get('database', '')} {plot.get('plot_type', '')}",
                    Path(plot["png"]),
                ))
    rseqc_images = [
        ("RSeQC gene-body coverage (5-prime to 3-prime)",
         results / "01_qc" / "rseqc" / "gene_body_coverage.svg"),
    ]
    sample_fields = [field for field in (
        "sample_id", "description", "genome", "condition", "batch", "subject",
        "biological_replicate_id", "technical_replicate_count", "sequencing_lane_count",
    ) if any(field in row for row in samples)]

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
        "", "## RSeQC annotation-aware QC", "",
        "RSeQC uses the exact selected annotation converted to BED12 (or the configured assembly-matched BED12).",
        "For QuantSeq REV, enrichment toward the transcript 3-prime end is expected; this is assay behavior, not conventional whole-transcript uniformity.",
        "",
        "| Sample | Genome | Condition | Dominant orientation | 5-prime fraction | 3-prime fraction | 3:5 ratio | 3-prime enriched |",
        "|---|---|---|---|---:|---:|---:|---|",
    ]
    lines.extend(
        f"| {row.get('sample_id', '')} | {row.get('genome', '')} | {row.get('condition', '')} | "
        f"{row.get('dominant_orientation', '')} | {row.get('five_prime_fraction', '')} | "
        f"{row.get('three_prime_fraction', '')} | {row.get('three_to_five_ratio', '')} | "
        f"{row.get('quantseq_three_prime_enriched', '')} |" for row in rseqc_rows
    )
    if not rseqc_rows:
        lines.append("| unavailable |  |  |  |  |  |  |  |")
    lines += [
        "", "## FastQ Screen contamination/species QC", "",
        f"- Status: `{fastq_screen_status}`",
        "- Per-lane and per-mate inventory: `../01_qc/fastq_screen/fastq_screen_summary.tsv`",
    ]
    lines += [
        "", "## Differential and APA summary", "",
        "| Contrast | DGE significant | DGE up | DGE down | APA-A significant genes | APA-A significant sites | APA-B significant genes | APA-B confirmed sites | PCPA A/B | Agreement % |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    lines.extend(
        f"| {row['contrast_id']} | {row['dge_significant']} | {row['dge_up']} | {row['dge_down']} | "
        f"{row['apa_a_significant_genes']} | {row['apa_a_significant_sites']} | "
        f"{row['apa_b_significant_genes']} | {row['apa_b_confirmed_sites']} | "
        f"{row['apa_a_pcpa']}/{row['apa_b_pcpa']} | {row['apa_direction_agreement_pct']} |"
        for row in contrast_rows
    )
    lines += [
        "", "Detailed, sortable result summaries:", "",
        "- `differential_gene_expression_summary.tsv` and `top_differential_genes.tsv`",
        "- `alternative_polyadenylation_summary.tsv` and `top_apa_gene_events.tsv`",
        "- `top_enrichment_terms.tsv`",
    ]
    lines += ["", "## Gene-set enrichment", "",
              "| Analysis | Genome | Contrast | Foreground genes | Significant ORA | Significant GSEA | Rich plots |",
              "|---|---|---|---:|---:|---:|---:|"]
    lines.extend(
        f"| {row.get('analysis_type', '')} | {row.get('genome', '')} | {row.get('contrast_id', '')} | "
        f"{row.get('foreground_genes', '')} | {row.get('significant_ora_terms', '')} | "
        f"{row.get('significant_gsea_terms', '')} | {row.get('rich_plot_count', '')} |" for row in enrichment_rows
    )
    if not enrichment_rows:
        lines.append("| disabled/unavailable |  |  |  |  |  |  |")
    lines += ["", "## APA-B validation and interpretation", "", f"- Status: `{apa_b_status}`"]
    lines.extend(f"- {row['property']}: {row['value']}" for row in apa_b_rows)
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
    lines += ["", "## Track collections", "", "| Collection | Description | BigWigs |", "|---|---|---:|"]
    lines.extend(
        f"| {row['collection']} | {row['description']} | {row['bigwigs']} |" for row in track_rows
    )
    lines += [
        "", "## Main files", "",
        "- `bigwig_collections.txt`: one-column, collection-grouped BigWig inventory.",
        "- `ucsc_track_descriptors/UCSC_bigWig_tracks.oneline.txt`: combined UCSC custom-track lines.",
        "- `ucsc_track_descriptors/<family>.txt`: one-line descriptors separated by normalization group.",
        "- `ucsc_track_descriptors/UCSC_trackDb.txt`: compatibility copy using the same one-line syntax.",
        "- `ucsc_track_descriptors/UCSC_descriptor_validation.tsv`: syntax-validation audit.",
        "- `IGV_session.xml`: session containing every generated BigWig.",
        "- `provenance_dashboard/`: receipt, environment, reference and complete output inventories.",
    ]
    markdown.write_text("\n".join(lines) + "\n", encoding="utf-8")

    c_legend = [
        {"stage": "C0", "meaning": "eligible end-defining molecules (R1 fragments in PE)", "use": "exact-end funnel; BAM/all-read tracks retain both PE mates"},
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
        {"property": "Generated BigWigs", "value": sum(int(row["bigwigs"]) for row in track_rows)},
        {"property": "Track collections", "value": len(track_rows)},
    ]
    artifact_html = "<ul>" + "".join(
        f"<li><a href='{html.escape(row['path'])}'>{html.escape(row['artifact'])}</a></li>"
        for row in artifact_rows
    ) + "</ul>"
    body = "".join([
        "<nav><a href='#overview'>Overview</a> | <a href='#samples'>Samples</a> | <a href='#qc'>QC</a> | "
        "<a href='#results'>Differential results</a> | <a href='#plots'>Plots</a> | "
        "<a href='#enrichment'>Enrichment</a> | <a href='#apa-b'>APA-B</a> | "
        "<a href='#tracks'>Tracks</a> | <a href='#provenance'>Provenance</a></nav>",
        f"<h1>rna_ends2tracks report: {html.escape(plan.project['project_id'])}</h1>",
        "<p class='note'>QuantSeq REV gene-expression and alternative-polyadenylation summary. "
        "Intragenic PAS are candidate premature cleavage/polyadenylation events; this assay alone "
        "does not prove loss of downstream engaged RNA polymerase.</p>",
        "<h2 id='overview'>Run overview</h2>", _html_table(overview, ["property", "value"]),
        "<h2>Count-universe legend</h2>", _html_table(c_legend, ["stage", "meaning", "use"]),
        "<h2>Module status</h2>", _html_table(module_rows, ["module", "status"]),
        "<h2 id='samples'>Validated samples</h2>",
        _html_table(samples, sample_fields) if sample_fields else "<p>No sample table available.</p>",
        "<h2 id='qc'>Sequencing, alignment, and 3-prime-end QC</h2>",
        "<p>The STAR table reports each technical-library/lane alignment. The orientation table verifies "
        "the reverse-compatible QuantSeq REV signal. In PE, both mate alignments remain in the BAM but the "
        "end-analysis C0 funnel counts R1 fragments. Count-universe funnels must satisfy "
        "<code>C0 = C1 + C1S</code> and <code>C1 = C2 + C2R</code>.</p>",
        "<h3>FastQ Screen contamination/species check</h3>",
        "<p>FastQ Screen examines a configured subset of every lane and mate against the site's species and "
        "contaminant databases. A skipped status means that no readable FASTQ_SCREEN_CONFIG was supplied; it "
        "does not mean that contamination was absent. PASS means the screen completed technically; database "
        "percentages still require interpretation for the expected species and local contaminant panel.</p>",
        _html_table(fastq_screen_rows, ["sample_id", "technical_replicate_id", "lane_id", "layout",
                                        "mates", "status", "text_reports", "config"])
        if fastq_screen_rows else "<p>FastQ Screen information was unavailable.</p>",
        "<h4>Per-database mapping percentages</h4>",
        _html_table(fastq_screen_metrics, ["sample_id", "technical_replicate_id", "lane_id", "mate",
                                           "database", "reads_processed", "pct_unmapped",
                                           "pct_one_hit_one_library", "pct_multiple_hits_one_library",
                                           "pct_one_hit_multiple_libraries"])
        if fastq_screen_metrics else "<p>No per-database FastQ Screen metrics were available.</p>",
        "<h3>STAR mapping</h3>",
        _html_table(star_rows, ["sample_id", "lane", "input_reads", "uniquely_mapped_reads",
                                "uniquely_mapped_pct", "multimapped_pct", "unmapped_too_short_pct"])
        if star_rows else "<p>STAR final logs were not available.</p>",
        "<h3>QuantSeq REV orientation</h3>",
        _html_table(orientation_rows, [field for field in (
            "sample_id", "technical_replicate_id", "lane_id", "star_forward_count",
            "star_reverse_count", "reverse_compatible_fraction", "status",
        ) if any(field in row for row in orientation_rows)])
        if orientation_rows else "<p>Orientation audit was not available.</p>",
        "<h3>RSeQC annotation-aware QC</h3>",
        "<p><code>infer_experiment.py</code> independently checks library orientation; "
        "<code>read_distribution.py</code> reports annotated genomic-feature distribution; and "
        "<code>geneBody_coverage.py</code> measures relative coverage from transcript 5-prime to 3-prime. "
        "QuantSeq REV should be strongly 3-prime enriched, so a non-uniform gene-body profile is expected.</p>",
        _html_table(rseqc_rows, [field for field in (
            "sample_id", "genome", "condition", "infer_undetermined_fraction",
            "infer_forward_fraction", "infer_reverse_fraction", "dominant_orientation",
            "five_prime_fraction", "three_prime_fraction", "three_to_five_ratio",
            "quantseq_three_prime_enriched",
        ) if any(field in row for row in rseqc_rows)])
        if rseqc_rows else "<p>RSeQC was disabled or its outputs were unavailable.</p>",
        _image_gallery(outdir, rseqc_images),
        "<h3>Exact-end filtering funnel</h3>",
        _html_table(funnel_rows, ["sample_id", "genome", "C0", "C1", "C1S", "C2", "C2R",
                                  "C1_over_C0_pct", "C2_over_C1_pct", "mask_rescued"])
        if funnel_rows else "<p>Exact-end audits were not available.</p>",
        "<h3>Active-PAS assignment</h3>",
        _html_table(active_pas_rows, ["genome", "samples", "active_pas", "C2_total", "C3_total",
                                     "C3_over_C2_pct", "ambiguous_pas", "pcpa_candidate_sites"])
        if active_pas_rows else "<p>Active-PAS audits were not available.</p>",
        "<h2 id='results'>Contrast-level results</h2>",
        "<p>Blank cells mean that the corresponding optional analysis was disabled or unavailable. "
        "Counts are recalculated from the result tables and checked against each module index before the report is published.</p>",
        "<label>Filter contrasts: <input id='contrast-filter' type='search' placeholder='condition, genome, result...'></label>",
        _html_table(contrast_rows, CONTRAST_SUMMARY_FIELDS, "contrast-table"),
        "<h3>Top differential genes</h3>",
        "<p>Up to 25 FDR-significant genes per contrast, ordered by adjusted p-value and absolute fold change. "
        "The complete DESeq2 tables remain linked in the output manifest.</p>",
        _html_table(dge_events, DGE_EVENT_FIELDS) if dge_events else "<p>No DGE gene passed the configured FDR.</p>",
        "<h3>Top APA gene-level events</h3>",
        "<p>Up to 25 FDR-significant genes per method and contrast. APA-A and APA-B remain independent; "
        "candidate PCPA marks intragenic premature cleavage/polyadenylation candidates, not proven termination.</p>",
        _html_table(apa_events, APA_EVENT_FIELDS) if apa_events else "<p>No APA gene passed the configured FDR.</p>",
        "<h2 id='plots'>DGE exploratory and contrast plots</h2>",
        "<p>PCA and sample-distance plots use variance-stabilized C4 counts. MA and volcano plots are "
        "generated independently for every pairwise contrast.</p>",
        _image_gallery(outdir, plot_images),
        "<h2 id='enrichment'>Gene-set enrichment</h2>",
        "<p>ORA uses significant foreground genes; ranked GSEA uses the complete tested background. "
        "APA queries distinguish any APA, distal/proximal shifts, and candidate PCPA increases/decreases. "
        "GO, Reactome, Hallmark, and KEGG collections are analyzed separately. Database-specific dotplots, "
        "barplots, and concept networks are indexed for every analysis. Mouse gene sets use the orthology "
        "mapping recorded in each provenance file.</p>",
        _html_table(enrichment_rows, [field for field in (
            "analysis_type", "genome", "contrast_id", "background_genes", "foreground_genes",
            "significant_ora_terms", "significant_gsea_terms", "rich_plot_count", "status",
        ) if any(field in row for row in enrichment_rows)]) if enrichment_rows else "<p>Enrichment was disabled or unavailable.</p>",
        "<h3>Top significant enrichment terms</h3>",
        _html_table(enrichment_terms, ENRICHMENT_TERM_FIELDS)
        if enrichment_terms else "<p>No enrichment term passed the configured adjusted-p threshold.</p>",
        _image_gallery(outdir, enrichment_images),
        "<h2 id='apa-b'>APA-B validation and interpretation</h2>",
        f"<p><strong>{html.escape(apa_b_status)}</strong></p>", _html_table(apa_b_rows, ["property", "value"]),
        "<p>APA-A and APA-B are independent analyses. Their catalogs are never merged; when both are validated, "
        "the workflow reports proximity and effect-direction concordance separately.</p>",
        "<h3>Top validated APA-B gene-level events</h3>",
        _html_table(apa_b_gene_events, [
            "genome", "contrast_id", "gene_id", "gene_padj", "shift", "max_abs_delta_PAU",
            "weighted_transcript_position_shift_nt", "confirmed_sites",
        ]) if apa_b_gene_events else "<p>No validated APA-B gene-level event passed the configured FDR, or APA-B was disabled.</p>",
        "<h2>Warnings</h2>", _html_table(warning_rows, ["warning_code", "message"]),
        "<h2 id='tracks'>BigWig track collections</h2>",
        "<p>Every row is a folder collection. Plus and minus are separate transcript-strand tracks. "
        "Minus BigWigs contain negative values; the generated UCSC descriptors can display their magnitude "
        "above zero with <code>negateValues=on</code>.</p>",
        _html_table(track_rows, ["collection", "description", "bigwigs", "transcript_plus", "transcript_minus"]),
        "<h2 id='provenance'>Provenance dashboard and main files</h2>",
        "<p>The dashboard inventories stage receipts, software/environment packages, reference and PAS-atlas "
        "identities, configuration checksums, and every output file. Large files are validated by size and "
        "modification time; small files also receive SHA-256 checksums.</p>", artifact_html,
        "<h2>Report scope</h2><p>This report integrates run status, samples, count-universe definitions, "
        "validated DGE/APA contrast counts, embedded statistical plots, enrichment results, warnings, browser-track "
        "inventories, RSeQC orientation/distribution/gene-body QC, and complete provenance. MultiQC provides "
        "the detailed FastQC, alignment, and RSeQC source reports. "
        "APA-B interpretation is explicitly marked unavailable unless its pinned pilot validation passes.</p>",
    ])
    html_target.write_text(
        "<!doctype html><html><head><meta charset='utf-8'><title>rna_ends2tracks report</title>"
        "<style>body{font-family:system-ui,sans-serif;max-width:1500px;margin:2rem auto;padding:0 1rem;line-height:1.45}"
        "h1,h2{color:#17324d}.note{background:#eef6fb;border-left:4px solid #2673a6;padding:.8rem 1rem}"
        ".table-wrap{overflow-x:auto;margin-bottom:1.5rem}table{border-collapse:collapse;width:100%;font-size:.9rem}"
        "th,td{border:1px solid #ccd5dd;padding:.4rem .55rem;text-align:left;white-space:nowrap}"
        "th{background:#e9eff4;position:sticky;top:0}tbody tr:nth-child(even){background:#f8fafb}"
        "a{color:#145d91}nav{position:sticky;top:0;background:#fff;padding:.7rem 0;border-bottom:1px solid #ccd5dd}"
        "input[type=search]{min-width:22rem;padding:.35rem}.gallery{display:grid;grid-template-columns:repeat(auto-fit,minmax(310px,1fr));gap:1rem}"
        "figure{margin:0;border:1px solid #ccd5dd;padding:.6rem;background:#fff}figure img{width:100%;height:auto}figcaption{font-size:.85rem;margin-top:.4rem}</style></head>"
        f"<body>{body}<script>const q=document.getElementById('contrast-filter');"
        "if(q){q.addEventListener('input',()=>{const v=q.value.toLowerCase();"
        "document.querySelectorAll('#contrast-table tbody tr').forEach(r=>r.hidden=!r.textContent.toLowerCase().includes(v));});}"
        "</script></body></html>\n", encoding="utf-8")
    provenance_outputs = generate_provenance_dashboard(plan, results, outdir)
    outputs = [
        markdown, html_target, summary, contrast_summary, dge_summary, apa_summary,
        dge_events_path, apa_events_path, enrichment_terms_path,
        *browser, *provenance_outputs,
    ]
    write_receipt("report", outdir, signature, outputs, ["rna-ends2tracks", "report"])
    event(log_dir, "report", "completed", str(html_target))
    return html_target
