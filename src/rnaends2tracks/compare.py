from __future__ import annotations

import csv
from collections import Counter, defaultdict
from pathlib import Path

from .config import RunPlan, signature_for
from .external import event
from .receipts import receipt_valid, write_receipt


def _read(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _write(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)


def _crosswalk(a: list[dict[str, str]], b: list[dict[str, str]], tolerance: int) -> list[dict[str, object]]:
    index: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in b:
        index[(row["chrom"], row["strand"])].append(row)
    for rows in index.values():
        rows.sort(key=lambda row: (int(row["start"]), row["pas_id"]))
    used: set[str] = set()
    result: list[dict[str, object]] = []
    for left in a:
        candidates = [row for row in index[(left["chrom"], left["strand"])]
                      if row["pas_id"] not in used and abs(int(row["start"]) - int(left["start"])) <= tolerance]
        same_gene = [row for row in candidates if left.get("gene_id") and row.get("gene_id") == left.get("gene_id")]
        candidates = same_gene or candidates
        right = min(candidates, key=lambda row: (abs(int(row["start"]) - int(left["start"])), row["pas_id"])) if candidates else None
        if right:
            used.add(right["pas_id"])
        result.append({
            "apa_a_pas_id": left["pas_id"], "apa_b_pas_id": right["pas_id"] if right else "",
            "chrom": left["chrom"], "strand": left["strand"],
            "distance_nt": abs(int(right["start"]) - int(left["start"])) if right else "",
            "apa_a_gene_id": left.get("gene_id", ""), "apa_b_gene_id": right.get("gene_id", "") if right else "",
            "match_class": "matched" if right else "apa_a_only",
        })
    for right in b:
        if right["pas_id"] not in used:
            result.append({"apa_a_pas_id": "", "apa_b_pas_id": right["pas_id"], "chrom": right["chrom"],
                           "strand": right["strand"], "distance_nt": "", "apa_a_gene_id": "",
                           "apa_b_gene_id": right.get("gene_id", ""), "match_class": "apa_b_only"})
    return result


def _effects(genome: str, matches: list[dict[str, object]], a_index: Path, b_index: Path) -> list[dict[str, object]]:
    a_files = {row["contrast_id"]: Path(row["result_file"]) for row in _read(a_index)}
    b_files = {row["contrast_id"]: Path(row["result_file"]) for row in _read(b_index)}
    rows: list[dict[str, object]] = []
    for contrast in sorted(set(a_files) & set(b_files)):
        a = {row["pas_id"]: row for row in _read(a_files[contrast])}
        b = {row["feature_id"]: row for row in _read(b_files[contrast])}
        for match in matches:
            aid, bid = str(match["apa_a_pas_id"]), str(match["apa_b_pas_id"])
            if aid not in a or bid not in b:
                continue
            try:
                da, db = float(a[aid]["delta_PAU"]), float(b[bid]["delta_PAU"])
            except (KeyError, ValueError):
                continue
            rows.append({"genome": genome, "contrast_id": contrast, "apa_a_pas_id": aid,
                         "apa_b_pas_id": bid, "apa_a_delta_PAU": da, "apa_b_delta_PAU": db,
                         "direction_agreement": "zero_effect" if da == 0 or db == 0 else (da > 0) == (db > 0),
                         "distance_nt": match["distance_nt"]})
    return rows


def _pcpa_agreement(
    genome: str, matches: list[dict[str, object]], a_path: Path, b_path: Path,
) -> list[dict[str, object]]:
    a = {(row["contrast_id"], row["pas_id"]) for row in _read(a_path)}
    b = {(row["contrast_id"], row["pas_id"]) for row in _read(b_path)}
    contrasts = sorted({key[0] for key in a | b})
    rows: list[dict[str, object]] = []
    for match in matches:
        if match["match_class"] != "matched":
            continue
        aid, bid = str(match["apa_a_pas_id"]), str(match["apa_b_pas_id"])
        for contrast in contrasts:
            in_a, in_b = (contrast, aid) in a, (contrast, bid) in b
            if in_a or in_b:
                rows.append({"genome": genome, "contrast_id": contrast, "apa_a_pas_id": aid,
                             "apa_b_pas_id": bid, "apa_a_candidate": in_a,
                             "apa_b_candidate": in_b, "agreement": in_a and in_b})
    return rows


def compare_apa(plan: RunPlan, results: Path, tolerance: int = 24, force: bool = False) -> None:
    outdir = results / "08_apa_comparison"
    log_dir = results / "logs"
    inputs: list[Path] = []
    for genome in plan.references:
        inputs.extend([results / "04_active_pas" / genome / "active_pas_catalog.tsv",
                       results / "06_apa_a_mcell2019" / genome / "dexseq" / "result_index.tsv",
                       results / "06_apa_a_mcell2019" / genome / "candidate_pcpa.tsv",
                       results / "07_apa_b" / genome / "pas_catalog.tsv",
                       results / "07_apa_b" / genome / "drimseq" / "result_index.tsv",
                       results / "07_apa_b" / genome / "candidate_pcpa.tsv"])
    signature = signature_for(inputs, {"module": "apa_comparison", "tolerance": tolerance})
    if not force and receipt_valid(outdir, signature):
        event(log_dir, "apa_comparison", "skipped", "Valid matching receipt")
        return
    all_matches: list[dict[str, object]] = []
    all_effects: list[dict[str, object]] = []
    all_pcpa: list[dict[str, object]] = []
    for genome in plan.references:
        a_catalog = results / "04_active_pas" / genome / "active_pas_catalog.tsv"
        b_catalog = results / "07_apa_b" / genome / "pas_catalog.tsv"
        matches = _crosswalk(_read(a_catalog), _read(b_catalog), tolerance)
        for row in matches:
            row["genome"] = genome
        all_matches.extend(matches)
        all_effects.extend(_effects(genome, matches,
            results / "06_apa_a_mcell2019" / genome / "dexseq" / "result_index.tsv",
            results / "07_apa_b" / genome / "drimseq" / "result_index.tsv"))
        all_pcpa.extend(_pcpa_agreement(genome, matches,
            results / "06_apa_a_mcell2019" / genome / "candidate_pcpa.tsv",
            results / "07_apa_b" / genome / "candidate_pcpa.tsv"))
    crosswalk = outdir / "site_crosswalk.tsv"
    fields = ["genome", "apa_a_pas_id", "apa_b_pas_id", "chrom", "strand", "distance_nt",
              "apa_a_gene_id", "apa_b_gene_id", "match_class"]
    _write(crosswalk, all_matches, fields)
    effect_path = outdir / "effect_concordance.tsv"
    _write(effect_path, all_effects, ["genome", "contrast_id", "apa_a_pas_id", "apa_b_pas_id",
        "apa_a_delta_PAU", "apa_b_delta_PAU", "direction_agreement", "distance_nt"])
    counts = Counter((str(row["genome"]), str(row["match_class"])) for row in all_matches)
    summary_rows = [{"genome": genome, "class": cls, "site_count": count}
                    for (genome, cls), count in sorted(counts.items())]
    summary = outdir / "summary.tsv"
    _write(summary, summary_rows, ["genome", "class", "site_count"])
    pcpa_path = outdir / "pcpa_agreement.tsv"
    _write(pcpa_path, all_pcpa, ["genome", "contrast_id", "apa_a_pas_id", "apa_b_pas_id",
                                 "apa_a_candidate", "apa_b_candidate", "agreement"])
    write_receipt("apa_comparison", outdir, signature, [crosswalk, effect_path, summary, pcpa_path],
                  ["rna-ends2tracks", "apa_comparison"])
    event(log_dir, "apa_comparison", "completed", f"Compared independent catalogs within {tolerance} nt")
