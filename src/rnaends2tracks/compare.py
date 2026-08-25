from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path

from .config import signature_for
from .external import event
from .receipts import receipt_valid, write_receipt


def _read(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def compare_apa(results: Path, tolerance: int = 24, force: bool = False) -> None:
    outdir = results / "06_apa_comparison"
    log_dir = results / "provenance" / "logs"
    outdir.mkdir(parents=True, exist_ok=True)
    input_paths = [results / "04_apa_a_repository" / "pas_catalog.tsv",
                   results / "05_apa_b_polyaseqtrap_drimseq" / "pas_catalog.tsv",
                   results / "04_apa_a_repository" / "dexseq" / "result_index.tsv",
                   results / "05_apa_b_polyaseqtrap_drimseq" / "drimseq" / "result_index.tsv",
                   results / "04_apa_a_repository" / "candidate_pcpa.tsv",
                   results / "05_apa_b_polyaseqtrap_drimseq" / "candidate_pcpa.tsv"]
    for index_path in input_paths[2:4]:
        input_paths.extend(Path(row["result_file"]) for row in _read(index_path))
    signature = signature_for(input_paths, {"module": "compare", "tolerance": tolerance})
    if not force and receipt_valid(outdir, signature):
        event(log_dir, "compare", "skipped", "Valid matching receipt")
        return
    a = _read(input_paths[0])
    b = _read(input_paths[1])
    b_index: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in b:
        b_index[(row["chrom"], row["strand"])].append(row)
    for rows in b_index.values():
        rows.sort(key=lambda row: int(row["start"]))
    crosswalk: list[dict[str, object]] = []
    used_b: set[str] = set()
    for row_a in a:
        candidates = [row for row in b_index[(row_a["chrom"], row_a["strand"])]
                      if row["pas_id"] not in used_b and abs(int(row["start"]) - int(row_a["start"])) <= tolerance]
        if row_a.get("gene_id"):
            gene_candidates = [row for row in candidates if row.get("gene_id") == row_a["gene_id"]]
            candidates = gene_candidates or candidates
        if candidates:
            row_b = min(candidates, key=lambda row: (abs(int(row["start"]) - int(row_a["start"])), row["pas_id"]))
            used_b.add(row_b["pas_id"])
            crosswalk.append({
                "apa_a_pas_id": row_a["pas_id"], "apa_b_pas_id": row_b["pas_id"],
                "chrom": row_a["chrom"], "strand": row_a["strand"],
                "distance_nt": abs(int(row_b["start"]) - int(row_a["start"])),
                "apa_a_gene_id": row_a.get("gene_id", ""), "apa_b_gene_id": row_b.get("gene_id", ""),
                "match_class": "matched",
            })
        else:
            crosswalk.append({
                "apa_a_pas_id": row_a["pas_id"], "apa_b_pas_id": "", "chrom": row_a["chrom"],
                "strand": row_a["strand"], "distance_nt": "", "apa_a_gene_id": row_a.get("gene_id", ""),
                "apa_b_gene_id": "", "match_class": "apa_a_only",
            })
    for row_b in b:
        if row_b["pas_id"] not in used_b:
            crosswalk.append({
                "apa_a_pas_id": "", "apa_b_pas_id": row_b["pas_id"], "chrom": row_b["chrom"],
                "strand": row_b["strand"], "distance_nt": "", "apa_a_gene_id": "",
                "apa_b_gene_id": row_b.get("gene_id", ""), "match_class": "apa_b_only",
            })
    path = outdir / "site_crosswalk.tsv"
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(crosswalk[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(crosswalk)
    summary = defaultdict(int)
    for row in crosswalk:
        summary[str(row["match_class"])] += 1
    with (outdir / "summary.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["class", "site_count"]); writer.writerows(sorted(summary.items()))
    _compare_effects(results, outdir, crosswalk)
    _compare_pcpa(results, outdir, crosswalk)
    write_receipt("compare", outdir, signature,
                  [outdir / "site_crosswalk.tsv", outdir / "summary.tsv", outdir / "effect_concordance.tsv", outdir / "pcpa_agreement.tsv"],
                  ["rna-ends2tracks", "compare"])
    event(log_dir, "compare", "completed", f"Matched APA sites within {tolerance} nt without merging catalogs")


def _compare_effects(results: Path, outdir: Path, crosswalk: list[dict[str, object]]) -> None:
    index_a = {row["contrast_id"]: row["result_file"] for row in _read(results / "04_apa_a_repository" / "dexseq" / "result_index.tsv")}
    index_b = {row["contrast_id"]: row["result_file"] for row in _read(results / "05_apa_b_polyaseqtrap_drimseq" / "drimseq" / "result_index.tsv")}
    output: list[dict[str, object]] = []
    for contrast in sorted(set(index_a) & set(index_b)):
        a = {row["pas_id"]: row for row in _read(Path(index_a[contrast]))}
        b = {row["feature_id"]: row for row in _read(Path(index_b[contrast]))}
        for match in crosswalk:
            a_id, b_id = str(match["apa_a_pas_id"]), str(match["apa_b_pas_id"])
            if not a_id or not b_id or a_id not in a or b_id not in b:
                continue
            try:
                da = float(a[a_id].get("delta_PAU", "nan")); db = float(b[b_id].get("delta_PAU", "nan"))
            except ValueError:
                continue
            output.append({
                "contrast_id": contrast, "apa_a_pas_id": a_id, "apa_b_pas_id": b_id,
                "apa_a_delta_PAU": da, "apa_b_delta_PAU": db,
                "direction_agreement": (da > 0) == (db > 0) if da != 0 and db != 0 else "zero_effect",
                "distance_nt": match["distance_nt"],
            })
    headers = ["contrast_id", "apa_a_pas_id", "apa_b_pas_id", "apa_a_delta_PAU", "apa_b_delta_PAU", "direction_agreement", "distance_nt"]
    with (outdir / "effect_concordance.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=headers, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(output)


def _compare_pcpa(results: Path, outdir: Path, crosswalk: list[dict[str, object]]) -> None:
    pcpa_a = {(row["contrast_id"], row["pas_id"]): row for row in _read(results / "04_apa_a_repository" / "candidate_pcpa.tsv")}
    pcpa_b = {(row["contrast_id"], row["pas_id"]): row for row in _read(results / "05_apa_b_polyaseqtrap_drimseq" / "candidate_pcpa.tsv")}
    rows: list[dict[str, object]] = []
    for match in crosswalk:
        if match["match_class"] != "matched":
            continue
        for contrast in sorted({key[0] for key in pcpa_a} | {key[0] for key in pcpa_b}):
            in_a = (contrast, str(match["apa_a_pas_id"])) in pcpa_a
            in_b = (contrast, str(match["apa_b_pas_id"])) in pcpa_b
            if in_a or in_b:
                rows.append({"contrast_id": contrast, "apa_a_pas_id": match["apa_a_pas_id"],
                             "apa_b_pas_id": match["apa_b_pas_id"], "apa_a_candidate": in_a,
                             "apa_b_candidate": in_b, "agreement": in_a and in_b})
    headers = ["contrast_id", "apa_a_pas_id", "apa_b_pas_id", "apa_a_candidate", "apa_b_candidate", "agreement"]
    with (outdir / "pcpa_agreement.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=headers, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)
