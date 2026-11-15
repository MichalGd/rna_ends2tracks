"""Explicit expression method identities and output discovery."""
from __future__ import annotations

import csv
from pathlib import Path

NEW_METHOD = "extended_exons_v1"
LEGACY_METHOD = "legacy_c4"
MODE_ORDER = {"stringent": 0, "relaxed_fdr": 1, "shared_gsea": 2, "legacy_c4": 3}


def expression_methods(project: dict) -> list[str]:
    # Pre-alpha13 YAML/RunPlan objects retain their historical contract. The
    # public .conf parser explicitly supplies the new default.
    return list(project.get("gene_expression", {}).get("methods", [LEGACY_METHOD]))


def expression_index(results: Path, genome: str, method: str) -> Path:
    leaf = "C4_primary_deseq2" if method == LEGACY_METHOD else NEW_METHOD
    return results / "05_gene_expression" / genome / leaf / "result_index.tsv"


def expression_rows(project: dict, results: Path, genome: str) -> list[dict[str, str]]:
    rows = []
    for method in expression_methods(project):
        path = expression_index(results, genome, method)
        if not path.is_file():
            continue
        with path.open(encoding="utf-8-sig", newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                row.setdefault("expression_method", method)
                row.setdefault("analysis_mode", "legacy_c4")
                settings = project.get("gene_expression", {})
                if row["analysis_mode"] == "stringent" and not settings.get("stringent", True):
                    continue
                if row["analysis_mode"] == "relaxed_fdr" and not settings.get("relaxed", True):
                    continue
                row["genome"] = genome
                rows.append(row)
    return sorted(rows, key=lambda row: (MODE_ORDER.get(row["analysis_mode"], 9), row["contrast_id"]))


def sample_bam(project: dict, results: Path, sample_id: str) -> Path:
    root = Path(project["alignment_source_dir"]) if project.get("input_mode") == "bam" else results / "02_alignment"
    return root / sample_id / f"{sample_id}.bam"
