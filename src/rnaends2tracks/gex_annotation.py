"""Deterministic exon-union annotation with bounded outer 3-prime extension.

Internal coordinates are zero-based, half-open. SAF export is one-based inclusive.
The original GTF and all original exon features are preserved.
"""
from __future__ import annotations

import bisect
import csv
import gzip
import json
import re
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path

from .receipts import sha256

POLICY_VERSION = 1
EXTEND_BIOTYPES = frozenset({
    "protein_coding", "lncRNA", "lincRNA", "antisense", "processed_transcript",
    "sense_intronic", "sense_overlapping", "3prime_overlapping_ncRNA", "macro_lncRNA",
    "IG_C_gene", "IG_D_gene", "IG_J_gene", "IG_V_gene", "TR_C_gene", "TR_D_gene", "TR_J_gene", "TR_V_gene",
})


@dataclass
class Model:
    gene_id: str
    chrom: str
    strand: str
    start: int
    end: int
    name: str = ""
    biotype: str = ""
    exons: list[tuple[int, int]] = field(default_factory=list)
    tss: set[int] = field(default_factory=set)


def merge_intervals(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
    merged: list[tuple[int, int]] = []
    for start, end in sorted(intervals):
        if merged and start <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(end, merged[-1][1]))
        else:
            merged.append((start, end))
    return merged


class Intervals:
    def __init__(self, records):
        self.records = sorted(records)
        self.starts = [row[0] for row in self.records]
        self.max_ends = []
        largest = -1
        for row in self.records:
            largest = max(largest, row[1])
            self.max_ends.append(largest)

    def overlapping(self, start, end):
        left = bisect.bisect_right(self.max_ends, start)
        right = bisect.bisect_left(self.starts, end)
        return (row for row in self.records[left:right] if row[1] > start)


def load_models(gtf: Path, sizes: dict[str, int]) -> dict[str, Model]:
    models: dict[str, Model] = {}
    opener = gzip.open if gtf.suffix == ".gz" else open
    with opener(gtf, "rt", encoding="utf-8") as handle:
        for number, line in enumerate(handle, 1):
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                raise ValueError(f"Invalid GTF record at line {number}")
            chrom, _, feature, first, last, _, strand, _, attributes = fields
            if feature not in {"gene", "transcript", "exon"}:
                continue
            attrs = dict(re.findall(r'(\w+)\s+"([^\"]*)"', attributes))
            gene_id = attrs.get("gene_id")
            if not gene_id or strand not in {"+", "-"}:
                raise ValueError(f"Missing gene_id or strand at GTF line {number}")
            start, end = int(first) - 1, int(last)
            if chrom not in sizes or not 0 <= start < end <= sizes[chrom]:
                raise ValueError(f"GTF coordinates/contig incompatible with chromosome sizes at line {number}")
            model = models.setdefault(gene_id, Model(gene_id, chrom, strand, start, end))
            if (model.chrom, model.strand) != (chrom, strand):
                raise ValueError(f"Gene ID has conflicting contigs/strands: {gene_id}")
            model.start, model.end = min(model.start, start), max(model.end, end)
            model.name = attrs.get("gene_name", model.name)
            biotype = attrs.get("gene_type", attrs.get("gene_biotype", ""))
            if biotype:
                if model.biotype and model.biotype != biotype:
                    raise ValueError(f"Conflicting gene biotypes for {gene_id}")
                model.biotype = biotype
            if feature == "transcript":
                model.tss.add(start if strand == "+" else end - 1)
            if feature == "exon":
                model.exons.append((start, end))
    for model in models.values():
        model.exons = merge_intervals(model.exons)
    if not any(model.exons for model in models.values()):
        raise ValueError("GTF has no gene-associated exons")
    return models


def _clip(start, end, model, index):
    candidates = [row for row in index.overlapping(start, end) if row[2] != model.gene_id]
    if not candidates:
        return start, end, "", ""
    if model.strand == "+":
        barrier = min(candidates, key=lambda row: (max(start, row[0]), row[2], row[3]))
        end = max(start, barrier[0])
    else:
        barrier = min(candidates, key=lambda row: (-min(end, row[1]), row[2], row[3]))
        start = min(end, barrier[1])
    return start, end, barrier[3], barrier[2]


def build_annotation(gtf: Path, chrom_sizes: Path, outdir: Path, extension: int = 6000, buffer: int = 2000) -> dict[str, Path]:
    if extension < 0 or buffer < 0:
        raise ValueError("Extension and promoter buffer must be nonnegative")
    outdir.mkdir(parents=True, exist_ok=True)
    outputs = {key: outdir / name for key, name in {
        "saf": "extended_exons.saf", "bed": "extended_exons.bed",
        "genes": "gene_metadata.tsv", "audit": "extension_audit.tsv", "manifest": "annotation_manifest.json",
    }.items()}
    inputs = {"gtf_sha256": sha256(gtf), "chrom_sizes_sha256": sha256(chrom_sizes),
              "extension_max_nt": extension, "tss_buffer_nt": buffer, "policy_version": POLICY_VERSION,
              "extend_biotypes": sorted(EXTEND_BIOTYPES)}
    if outputs["manifest"].is_file():
        try:
            prior = json.loads(outputs["manifest"].read_text(encoding="utf-8"))
            if prior["inputs"] == inputs and set(prior["output_sha256"]) == set(outputs)-{"manifest"} and all(outputs[key].is_file() and sha256(outputs[key]) == value
                                                 for key, value in prior["output_sha256"].items()):
                return outputs
        except (OSError, ValueError, KeyError, TypeError):
            pass
    sizes = {}
    for line in chrom_sizes.read_text(encoding="utf-8").splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        chrom, size = line.split()[:2]
        if chrom in sizes or int(size) <= 0:
            raise ValueError(f"Invalid or duplicate chromosome size: {chrom}")
        sizes[chrom] = int(size)
    models = load_models(gtf, sizes)
    protected = defaultdict(list)
    for model in models.values():
        key = (model.chrom, model.strand)
        protected[key].append((model.start, model.end, model.gene_id, "other_gene_body"))
        for tss in model.tss or {model.start if model.strand == "+" else model.end - 1}:
            if buffer > 1:
                start, end = (max(0, tss - buffer + 1), tss) if model.strand == "+" else (tss + 1, min(sizes[model.chrom], tss + buffer))
                if start < end:
                    protected[key].append((start, end, model.gene_id, "tss_buffer"))
    barriers = {key: Intervals(rows) for key, rows in protected.items()}
    bodies = defaultdict(list)
    for model in models.values():
        bodies[(model.chrom, model.strand)].append((model.start, model.end, model.gene_id, "body"))
    body_indexes = {key: Intervals(rows) for key, rows in bodies.items()}
    candidates, audit = {}, {}
    for gene_id, model in sorted(models.items()):
        if not model.exons:
            continue
        anchor = model.exons[-1][1] if model.strand == "+" else model.exons[0][0]
        start, end = (anchor, min(sizes[model.chrom], anchor + extension)) if model.strand == "+" else (max(0, anchor - extension), anchor)
        reason = "full_cap" if end - start == extension else "chromosome_edge"
        neighbor = ""
        if not extension:
            reason = "extension_disabled"
        elif model.biotype not in EXTEND_BIOTYPES:
            start = end = anchor
            reason = "biotype_not_extended"
        elif start < end:
            start, end, limiting, neighbor = _clip(start, end, model, barriers[(model.chrom, model.strand)])
            reason = limiting or reason
        candidates[gene_id] = (start, end)
        audit[gene_id] = {"gene_id": gene_id, "chrom": model.chrom, "strand": model.strand,
                         "original_3prime": anchor if model.strand == "+" else anchor + 1,
                         "requested_extension_nt": extension, "reason": reason, "limiting_gene": neighbor,
                         "tss_source": "transcripts" if model.tss else "gene_boundary_fallback"}
    extension_groups = defaultdict(list)
    for gene_id, (start, end) in candidates.items():
        model = models[gene_id]
        if start < end:
            extension_groups[(model.chrom, model.strand)].append((start, end, gene_id, "extension_conflict"))
    extension_indexes = {key: Intervals(rows) for key, rows in extension_groups.items()}
    for gene_id, (start, end) in list(candidates.items()):
        model = models[gene_id]
        if start < end:
            start, end, reason, neighbor = _clip(start, end, model, extension_indexes[(model.chrom, model.strand)])
            if reason:
                audit[gene_id].update(reason=reason, limiting_gene=neighbor)
        candidates[gene_id] = start, end
        audit[gene_id].update(accepted_extension_nt=end-start,
                             extension_start=start+1 if start < end else "",
                             extension_end=end if start < end else "")
        opposite = body_indexes.get((model.chrom, "-" if model.strand == "+" else "+"))
        overlaps = sorted({row[2] for row in opposite.overlapping(start,end)}) if opposite and start < end else []
        audit[gene_id]["opposite_strand_extension_overlaps"] = ";".join(overlaps)
    with outputs["saf"].open("w", encoding="utf-8", newline="") as sf, outputs["bed"].open("w", encoding="utf-8", newline="") as bed, outputs["genes"].open("w", encoding="utf-8", newline="") as gf:
        writer = csv.writer(sf, delimiter="\t", lineterminator="\n")
        writer.writerow(["GeneID", "Chr", "Start", "End", "Strand"])
        genes = csv.writer(gf, delimiter="\t", lineterminator="\n")
        genes.writerow(["gene_id", "gene_name", "gene_biotype", "chrom", "strand", "original_exon_nt", "extended_exon_nt"])
        for gene_id in sorted(candidates):
            model = models[gene_id]
            start, end = candidates[gene_id]
            intervals = merge_intervals(model.exons + ([(start, end)] if start < end else []))
            for left, right in intervals:
                writer.writerow([gene_id, model.chrom, left+1, right, model.strand])
                bed.write(f"{model.chrom}\t{left}\t{right}\t{gene_id}\t0\t{model.strand}\n")
            genes.writerow([gene_id, model.name, model.biotype, model.chrom, model.strand,
                            sum(b-a for a,b in model.exons), sum(b-a for a,b in intervals)])
    with outputs["audit"].open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(next(iter(audit.values()))), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(audit[key] for key in sorted(audit))
    manifest = {"inputs": inputs, "gene_count": len(candidates),
                "output_sha256": {key: sha256(path) for key, path in outputs.items() if key != "manifest"},
                "coordinate_convention": "SAF: 1-based inclusive; BED: 0-based half-open"}
    outputs["manifest"].write_text(json.dumps(manifest, indent=2, sort_keys=True)+"\n", encoding="utf-8")
    return outputs
