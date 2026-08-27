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
    clean_description = f"{description}: {path.name}".replace('"', "'")
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
    combined = descriptor_dir / "UCSC_bigWig_tracks.oneline.txt"
    generated: list[Path] = [inventory, combined]
    all_lines: list[str] = []
    track_index = 0
    track_numbers: dict[Path, int] = {}
    for group_index, row in enumerate(collection_rows):
        group = str(row["collection"])
        description = str(row["description"])
        color = TRACK_COLLECTION_COLORS[group_index % len(TRACK_COLLECTION_COLORS)]
        group_lines = [f"# {group} - {description}"]
        for path in bigwigs:
            if path.relative_to(track_root).parent.as_posix() != group:
                continue
            track_index += 1
            track_numbers[path] = track_index
            group_lines.append(_ucsc_track_line(
                path, track_root, track_index, color, description, public_prefix,
                "../../09_tracks", negate_minus, view_limits,
            ))
        group_file = descriptor_dir / f"{group.replace('/', '__')}.txt"
        group_file.write_text("\n".join(group_lines) + "\n", encoding="utf-8")
        generated.append(group_file)
        all_lines.extend([*group_lines, ""])
    combined.write_text("\n".join(all_lines).rstrip() + "\n", encoding="utf-8")

    trackdb = outdir / "UCSC_trackDb.txt"
    legacy_lines: list[str] = []
    for group_index, row in enumerate(collection_rows):
        group = str(row["collection"])
        description = str(row["description"])
        color = TRACK_COLLECTION_COLORS[group_index % len(TRACK_COLLECTION_COLORS)]
        legacy_lines.extend([f"# {group} - {description}"])
        for path in bigwigs:
            if path.relative_to(track_root).parent.as_posix() != group:
                continue
            legacy_lines.append(_ucsc_track_line(
                path, track_root, track_numbers[path], color, description, public_prefix,
                "../09_tracks", negate_minus, view_limits,
            ))
        legacy_lines.append("")
    trackdb.write_text("\n".join(legacy_lines).rstrip() + "\n", encoding="utf-8")
    igv = outdir / "IGV_session.xml"
    resources = "\n".join(
        f'    <Resource path="../{html.escape(path.relative_to(results).as_posix())}"/>' for path in bigwigs
    )
    igv.write_text(f'<?xml version="1.0" encoding="UTF-8"?>\n<Session version="8">\n  <Resources>\n{resources}\n  </Resources>\n</Session>\n',
                   encoding="utf-8")
    generated.extend([trackdb, igv])
    return generated, collection_rows


def _main_artifacts(results: Path, outdir: Path) -> list[dict[str, str]]:
    candidates: list[tuple[str, Path]] = [
        ("MultiQC report", results / "01_qc" / "multiqc" / "multiqc_report.html"),
        ("Protocol-orientation audit", results / "02_alignment" / "protocol_orientation.tsv"),
        ("Track normalization table", results / "09_tracks" / "track_normalization.tsv"),
        ("One-column BigWig inventory", outdir / "bigwig_collections.txt"),
        ("Combined one-line UCSC descriptors", outdir / "ucsc_track_descriptors" / "UCSC_bigWig_tracks.oneline.txt"),
        ("IGV session", outdir / "IGV_session.xml"),
        ("Full contrast summary", outdir / "contrast_summary.tsv"),
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
    inputs.extend(sorted((results / "02_alignment").rglob("*Log.final.out")))
    inputs.append(results / "02_alignment" / "protocol_orientation.tsv")
    inputs.extend(sorted((results / "03_exact_ends").glob("*/*/end_audit.json")))
    inputs.extend(sorted((results / "04_active_pas").glob("*/count_universe_audit.json")))
    inputs = [path for path in inputs if path.is_file()]
    signature = signature_for(inputs, {
        "module": "report", "project": plan.project["project_id"],
        "reporting": plan.project.get("reporting", {}),
        "browser_tracks": {
            key: plan.project.get("tracks", {}).get(key)
            for key in ("ucsc_bigdata_url_prefix", "ucsc_negate_minus_tracks", "ucsc_view_limits")
        },
    })
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
    browser, track_rows = _browser_assets(plan, results, outdir)
    artifact_rows = _main_artifacts(results, outdir)
    samples = _rows(results / "00_metadata" / "validated_samples.tsv")
    star_rows = _star_qc_rows(results)
    orientation_rows = _rows(results / "02_alignment" / "protocol_orientation.tsv")
    funnel_rows = _exact_funnel_rows(results)
    active_pas_rows = _active_pas_rows(results)
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
    lines += ["", "## Track collections", "", "| Collection | Description | BigWigs |", "|---|---|---:|"]
    lines.extend(
        f"| {row['collection']} | {row['description']} | {row['bigwigs']} |" for row in track_rows
    )
    lines += [
        "", "## Main files", "",
        "- `bigwig_collections.txt`: one-column, collection-grouped BigWig inventory.",
        "- `ucsc_track_descriptors/UCSC_bigWig_tracks.oneline.txt`: combined UCSC custom-track lines.",
        "- `UCSC_trackDb.txt`: compatibility copy using the same one-line syntax.",
        "- `IGV_session.xml`: session containing every generated BigWig.",
    ]
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
        {"property": "Generated BigWigs", "value": sum(int(row["bigwigs"]) for row in track_rows)},
        {"property": "Track collections", "value": len(track_rows)},
    ]
    artifact_html = "<ul>" + "".join(
        f"<li><a href='{html.escape(row['path'])}'>{html.escape(row['artifact'])}</a></li>"
        for row in artifact_rows
    ) + "</ul>"
    body = "".join([
        "<nav><a href='#overview'>Overview</a> | <a href='#samples'>Samples</a> | <a href='#qc'>QC</a> | "
        "<a href='#results'>Differential results</a> | <a href='#tracks'>Tracks</a> | "
        "<a href='#files'>Main files</a></nav>",
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
        "the reverse-compatible QuantSeq REV signal. Count-universe funnels must satisfy "
        "<code>C0 = C1 + C1S</code> and <code>C1 = C2 + C2R</code>.</p>",
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
        "<h3>Exact-end filtering funnel</h3>",
        _html_table(funnel_rows, ["sample_id", "genome", "C0", "C1", "C1S", "C2", "C2R",
                                  "C1_over_C0_pct", "C2_over_C1_pct", "mask_rescued"])
        if funnel_rows else "<p>Exact-end audits were not available.</p>",
        "<h3>Active-PAS assignment</h3>",
        _html_table(active_pas_rows, ["genome", "samples", "active_pas", "C2_total", "C3_total",
                                     "C3_over_C2_pct", "ambiguous_pas", "pcpa_candidate_sites"])
        if active_pas_rows else "<p>Active-PAS audits were not available.</p>",
        "<h2 id='results'>Contrast-level results</h2>",
        "<p>Blank cells mean that the corresponding optional analysis was disabled or unavailable.</p>",
        "<label>Filter contrasts: <input id='contrast-filter' type='search' placeholder='condition, genome, result...'></label>",
        _html_table(contrast_rows, CONTRAST_SUMMARY_FIELDS, "contrast-table"),
        "<h2>Warnings</h2>", _html_table(warning_rows, ["warning_code", "message"]),
        "<h2 id='tracks'>BigWig track collections</h2>",
        "<p>Every row is a folder collection. Plus and minus are separate transcript-strand tracks. "
        "Minus BigWigs contain negative values; the generated UCSC descriptors can display their magnitude "
        "above zero with <code>negateValues=on</code>.</p>",
        _html_table(track_rows, ["collection", "description", "bigwigs", "transcript_plus", "transcript_minus"]),
        "<h2 id='files'>Main reports and data indexes</h2>", artifact_html,
        "<h2>Report scope</h2><p>This alpha.10 report integrates run status, samples, count-universe definitions, "
        "validated DGE/APA contrast counts, warnings, browser-track inventories, and links to primary indexes. "
        "MultiQC remains the detailed sequencing/alignment QC report. Enrichment plots and the planned "
        "independent APA-B interpretation are not fabricated when those modules are unavailable.</p>",
    ])
    html_target.write_text(
        "<!doctype html><html><head><meta charset='utf-8'><title>rna_ends2tracks report</title>"
        "<style>body{font-family:system-ui,sans-serif;max-width:1500px;margin:2rem auto;padding:0 1rem;line-height:1.45}"
        "h1,h2{color:#17324d}.note{background:#eef6fb;border-left:4px solid #2673a6;padding:.8rem 1rem}"
        ".table-wrap{overflow-x:auto;margin-bottom:1.5rem}table{border-collapse:collapse;width:100%;font-size:.9rem}"
        "th,td{border:1px solid #ccd5dd;padding:.4rem .55rem;text-align:left;white-space:nowrap}"
        "th{background:#e9eff4;position:sticky;top:0}tbody tr:nth-child(even){background:#f8fafb}"
        "a{color:#145d91}nav{position:sticky;top:0;background:#fff;padding:.7rem 0;border-bottom:1px solid #ccd5dd}"
        "input[type=search]{min-width:22rem;padding:.35rem}</style></head>"
        f"<body>{body}<script>const q=document.getElementById('contrast-filter');"
        "if(q){q.addEventListener('input',()=>{const v=q.value.toLowerCase();"
        "document.querySelectorAll('#contrast-table tbody tr').forEach(r=>r.hidden=!r.textContent.toLowerCase().includes(v));});}"
        "</script></body></html>\n", encoding="utf-8")
    outputs = [markdown, html_target, summary, contrast_summary, *browser]
    write_receipt("report", outdir, signature, outputs, ["rna-ends2tracks", "report"])
    event(log_dir, "report", "completed", str(html_target))
    return html_target
