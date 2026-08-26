#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import subprocess
import sys
import tempfile
from collections import Counter
from collections.abc import Iterable
from dataclasses import dataclass
from datetime import datetime, timezone
from itertools import chain
from pathlib import Path

# Permit administrator use directly from a source checkout as documented, while
# installed copies continue to import the packaged module normally.
SOURCE_ROOT = Path(__file__).resolve().parents[1] / "src"
if SOURCE_ROOT.is_dir():
    sys.path.insert(0, str(SOURCE_ROOT))

from rnaends2tracks.mcell2019 import (
    Site,
    assign_gene,
    build_gene_bins,
    load_gene_models,
    load_rescue_sites,
    read_chrom_sizes,
    rescue_overlap,
)


@dataclass(frozen=True)
class Record:
    chrom: str
    position: int
    strand: str
    source: str
    source_id: str
    gene_id: str = ""
    internal_priming: bool = False


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def opener(path: Path):
    return gzip.open(path, "rt", encoding="utf-8") if path.suffix == ".gz" else path.open(encoding="utf-8")


def read_records(path: Path, source: str) -> list[Record]:
    """Read BED6 or headered TSV; coordinates must already be 0-based half-open."""
    records: list[Record] = []
    with opener(path) as handle:
        lines = (line for line in handle if line.strip() and not line.startswith(("#", "track", "browser")))
        first = next(lines, None)
        if first is None:
            return records
        fields = first.rstrip("\n").split("\t")
        headered = fields[0].lower() in {"chrom", "chromosome", "seqnames"}
        header = {name.lower(): index for index, name in enumerate(fields)} if headered else {}
        iterable: Iterable[list[str]] = (line.rstrip("\n").split("\t") for line in lines)
        if not headered:
            iterable = chain([fields], iterable)
        for number, row in enumerate(iterable, start=2 if headered else 1):
            try:
                if headered:
                    chrom = row[header.get("chrom", header.get("chromosome", header.get("seqnames", -1)))]
                    start = int(row[header["start"]])
                    end = int(row[header["end"]])
                    strand = row[header["strand"]]
                    identifier = row[header["pas_id"]] if "pas_id" in header else f"{source}_{number}"
                    gene = row[header["gene_id"]] if "gene_id" in header else ""
                    flag = row[header["internal_priming_flag"]].lower() if "internal_priming_flag" in header else "false"
                else:
                    if len(row) < 6:
                        raise ValueError("BED input needs at least six columns")
                    chrom, start, end, identifier, _, strand = row[:6]
                    start, end = int(start), int(end)
                    gene = row[6] if len(row) > 6 else ""
                    flag = row[7].lower() if len(row) > 7 else "false"
                if start < 0 or end <= start or strand not in {"+", "-"}:
                    raise ValueError("invalid coordinate or strand")
            except (KeyError, IndexError, ValueError) as exc:
                raise SystemExit(f"Invalid {source} record at {path}:{number}: {exc}") from exc
            position = end - 1 if strand == "+" else start
            records.append(Record(chrom, position, strand, source, identifier, gene,
                                  flag in {"1", "true", "yes", "likely", "internal_priming"}))
    return records


def lift_records(records: list[Record], chain: Path, lift_over: str) -> tuple[list[Record], Counter[str]]:
    audit: Counter[str] = Counter()
    with tempfile.TemporaryDirectory(prefix="rna_ends_pas_liftover_") as temporary:
        root = Path(temporary)
        source, mapped, unmapped = root / "source.bed", root / "mapped.bed", root / "unmapped.bed"
        with source.open("w", encoding="utf-8") as handle:
            for index, row in enumerate(records):
                handle.write(f"{row.chrom}\t{row.position}\t{row.position + 1}\tR{index}\t0\t{row.strand}\n")
        subprocess.run([lift_over, "-multiple", str(source), str(chain), str(mapped), str(unmapped)], check=True)
        hits: dict[int, list[tuple[str, int, str]]] = {}
        if mapped.is_file():
            with mapped.open(encoding="utf-8") as handle:
                for line in handle:
                    chrom, start, _, token, _, strand = line.rstrip("\n").split("\t")[:6]
                    hits.setdefault(int(token[1:]), []).append((chrom, int(start), strand))
        result: list[Record] = []
        for index, record in enumerate(records):
            mapped_hits = hits.get(index, [])
            if len(mapped_hits) != 1:
                audit["unmapped" if not mapped_hits else "multiple"] += 1
                continue
            chrom, position, strand = mapped_hits[0]
            if strand != record.strand:
                audit["strand_changed"] += 1
                continue
            result.append(Record(chrom, position, strand, record.source, record.source_id,
                                 record.gene_id, record.internal_priming))
            audit["unique_preserved"] += 1
    return result, audit


def cluster(records: list[Record], distance: int) -> list[list[Record]]:
    grouped: dict[tuple[str, str], list[Record]] = {}
    for record in records:
        grouped.setdefault((record.chrom, record.strand), []).append(record)
    clusters: list[list[Record]] = []
    for key in sorted(grouped):
        current: list[Record] = []
        prior = -10**30
        for record in sorted(grouped[key], key=lambda row: (row.position, row.source, row.source_id)):
            if current and record.position - prior > distance:
                clusters.append(current); current = []
            current.append(record); prior = record.position
        if current:
            clusters.append(current)
    return clusters


def write_gzip_tsv(path: Path, fields: list[str], rows: list[dict[str, object]]) -> None:
    with gzip.GzipFile(filename=str(path), mode="wb", mtime=0) as compressed:
        import io
        with io.TextIOWrapper(compressed, encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
            writer.writeheader(); writer.writerows(rows)


def write_gzip_bed(path: Path, rows: list[dict[str, object]]) -> None:
    with gzip.GzipFile(filename=str(path), mode="wb", mtime=0) as compressed:
        import io
        with io.TextIOWrapper(compressed, encoding="utf-8", newline="") as handle:
            for row in rows:
                handle.write("\t".join(str(row[key]) for key in
                    ("chrom", "start", "end", "pas_id", "score", "strand")) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser(description="Build a versioned, assembly-locked PAS mask-rescue atlas")
    parser.add_argument("--species", choices=("human", "mouse"), required=True)
    parser.add_argument("--assembly", choices=("GRCh38", "GRCm39"), required=True)
    parser.add_argument("--annotation-release", required=True)
    parser.add_argument("--atlas-id", required=True)
    parser.add_argument("--source-manifest", type=Path, required=True,
                        help="JSON source ledger with file,url,release,download_date,license,sha256")
    parser.add_argument("--gtf", type=Path, required=True)
    parser.add_argument("--chrom-sizes", type=Path, required=True)
    parser.add_argument("--gencode", type=Path, required=True, help="Assembly-native BED6/headered TSV")
    parser.add_argument("--polyadb-main", type=Path, required=True)
    parser.add_argument("--polyadb-max", type=Path)
    parser.add_argument("--polyasite", type=Path)
    parser.add_argument("--polyasite-assembly")
    parser.add_argument("--liftover-chain", type=Path)
    parser.add_argument("--lift-over", default="liftOver")
    parser.add_argument("--merge-distance", type=int, default=10)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if (args.species, args.assembly) not in {("human", "GRCh38"), ("mouse", "GRCm39")}:
        raise SystemExit("Species/assembly mismatch")
    paths = [args.gtf, args.chrom_sizes, args.gencode, args.polyadb_main, args.source_manifest]
    paths += [path for path in (args.polyadb_max, args.polyasite, args.liftover_chain) if path]
    for path in paths:
        if not path.is_file():
            raise SystemExit(f"Input does not exist: {path}")
    try:
        source_ledger = json.loads(args.source_manifest.read_text(encoding="utf-8"))
        source_entries = source_ledger["sources"]
    except (OSError, ValueError, KeyError, TypeError) as exc:
        raise SystemExit("--source-manifest must be JSON with a sources array") from exc
    source_paths = [args.gencode, args.polyadb_main] + [path for path in
        (args.polyadb_max, args.polyasite, args.liftover_chain) if path]
    by_file = {str(entry.get("file", "")): entry for entry in source_entries if isinstance(entry, dict)}
    required_metadata = ("url", "release", "download_date", "license", "sha256")
    for source_path in source_paths:
        entry = by_file.get(source_path.name)
        if entry is None or any(not str(entry.get(key, "")).strip() for key in required_metadata):
            raise SystemExit(f"Source manifest lacks complete metadata for {source_path.name}")
        if entry["sha256"].lower() != sha256(source_path):
            raise SystemExit(f"Source-manifest checksum mismatch: {source_path}")

    records = read_records(args.gencode, "GENCODE") + read_records(args.polyadb_main, "PolyA_DB_Main")
    if args.polyadb_max:
        records += read_records(args.polyadb_max, "PolyA_DB_Max")
    lift_audit: Counter[str] = Counter()
    if args.polyasite:
        if not args.polyasite_assembly:
            raise SystemExit("--polyasite requires an explicit --polyasite-assembly")
        polyasite = read_records(args.polyasite, "PolyASite")
        if args.polyasite_assembly and args.polyasite_assembly != args.assembly:
            if not args.liftover_chain:
                raise SystemExit("Foreign PolyASite assembly requires --liftover-chain")
            polyasite, lift_audit = lift_records(polyasite, args.liftover_chain, args.lift_over)
        records += polyasite

    lengths = read_chrom_sizes(args.chrom_sizes)
    rejected = Counter()
    valid: list[Record] = []
    for record in records:
        if record.internal_priming:
            rejected["internal_priming_flag"] += 1
        elif record.chrom not in lengths or not 0 <= record.position < lengths[record.chrom]:
            rejected["out_of_bounds_or_unknown_contig"] += 1
        else:
            valid.append(record)
    genes = load_gene_models(args.gtf)
    gene_bins = build_gene_bins(genes)
    master: list[dict[str, object]] = []
    for index, members in enumerate(cluster(valid, args.merge_distance), start=1):
        sources = sorted({row.source for row in members})
        preferred = [row.position for row in members if row.source == "PolyA_DB_Main"] or \
                    [row.position for row in members if row.source == "GENCODE"] or [row.position for row in members]
        position = min(preferred)
        chrom, strand = members[0].chrom, members[0].strand
        if "PolyA_DB_Main" in sources and "GENCODE" in sources:
            tier, confidence = "core", "A"
        elif "PolyA_DB_Main" in sources or "GENCODE" in sources:
            tier, confidence = "core", "B"
        else:
            tier, confidence = "rescue", "C"
        gene_id, feature, assignment = assign_gene(
            Site(chrom, strand, position, position, position + 1), genes, bins=gene_bins
        )
        master.append({"chrom": chrom, "start": position, "end": position + 1,
                       "pas_id": f"{args.atlas_id}_PAS_{index:08d}", "score": len(members), "strand": strand,
                       "species": args.species, "assembly": args.assembly, "tier": tier,
                       "confidence": confidence, "sources": ";".join(sources),
                       "source_ids": ";".join(sorted({row.source_id for row in members})),
                       "source_gene_ids": ";".join(sorted({row.gene_id for row in members if row.gene_id})),
                       "gene_id": gene_id, "feature_class": feature, "assignment_status": assignment})
    args.output.mkdir(parents=True, exist_ok=True)
    fields = ["chrom", "start", "end", "pas_id", "score", "strand", "species", "assembly", "tier",
              "confidence", "sources", "source_ids", "source_gene_ids", "gene_id", "feature_class", "assignment_status"]
    if not any(row["tier"] == "core" for row in master):
        raise SystemExit("Atlas has no core PAS after filtering")
    write_gzip_tsv(args.output / "master.tsv.gz", fields, master)
    write_gzip_bed(args.output / "core.bed.gz", [row for row in master if row["tier"] == "core"])
    write_gzip_bed(args.output / "rescue.bed.gz", [row for row in master if row["tier"] == "rescue"])
    installed_sites = load_rescue_sites(args.output, "core_plus_rescue")
    if any(not rescue_overlap(str(row["chrom"]), str(row["strand"]), int(row["start"]),
                              (installed_sites,), 20) for row in master):
        raise SystemExit("Synthetic mask-rescue self-test failed")
    provenance = {"atlas_id": args.atlas_id, "species": args.species, "assembly": args.assembly,
        "annotation_release": args.annotation_release, "coordinate_system": "0-based half-open",
        "built_at": datetime.now(timezone.utc).isoformat(), "merge_distance": args.merge_distance,
        "inputs": [{"path": str(path.resolve()), "sha256": sha256(path)} for path in paths],
        "source_ledger": source_ledger,
        "liftover": dict(lift_audit), "rejected": dict(rejected),
        "counts": dict(Counter(str(row["tier"]) for row in master)),
        "synthetic_mask_rescue_test": "PASS"}
    (args.output / "provenance.json").write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    report = args.output / "build_report.tsv"
    with report.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["metric", "value"])
        for key, value in sorted({**provenance["counts"], **dict(rejected), **{f"liftover_{k}": v for k, v in lift_audit.items()}}.items()):
            writer.writerow([key, value])
    (args.output / "SHA256SUMS").write_text("", encoding="utf-8")
    for filename in ("core.bed.gz", "rescue.bed.gz", "master.tsv.gz", "provenance.json", "build_report.tsv"):
        with (args.output / "SHA256SUMS").open("a", encoding="utf-8") as handle:
            handle.write(f"{sha256(args.output / filename)}  {filename}\n")


if __name__ == "__main__":
    main()
