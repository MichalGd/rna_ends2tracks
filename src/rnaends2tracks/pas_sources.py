from __future__ import annotations

import csv
import gzip
import hashlib
import json
import subprocess
import tempfile
import zipfile
from collections import Counter
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path


@dataclass(frozen=True)
class PasSourceRecord:
    chrom: str
    position: int
    strand: str
    source_id: str
    source_gene: str = ""


SOURCE_PROFILES = {
    "human": {
        "assembly": "GRCh38",
        "annotation": "GENCODE_v42",
        "polyadb_assembly": "hg38",
        "gencode_url": (
            "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/"
            "release_42/gencode.v42.polyAs.gtf.gz"
        ),
        "polyadb_url": "https://exon.apps.wistar.org/polya_db/v4/download/4.1/HumanPas.zip",
    },
    "mouse": {
        "assembly": "GRCm39",
        "annotation": "GENCODE_vM31",
        "polyadb_assembly": "mm10",
        "gencode_url": (
            "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/"
            "release_M31/gencode.vM31.polyAs.gtf.gz"
        ),
        "polyadb_url": "https://exon.apps.wistar.org/polya_db/v4/download/4.1/MousePas.zip",
    },
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _open_text(path: Path):
    return gzip.open(path, "rt", encoding="utf-8") if path.suffix == ".gz" else path.open(encoding="utf-8")


def _gtf_attributes(value: str) -> dict[str, str]:
    result: dict[str, str] = {}
    for field in value.split(";"):
        pieces = field.strip().split(" ", 1)
        if len(pieces) == 2:
            result[pieces[0]] = pieces[1].strip().strip('"')
    return result


def read_gencode_polya_sites(path: Path) -> tuple[list[PasSourceRecord], Counter[str]]:
    """Read only cleavage-site features; upstream signal motifs are not PAS coordinates."""
    records: dict[tuple[str, int, str], PasSourceRecord] = {}
    audit: Counter[str] = Counter()
    with _open_text(path) as handle:
        for number, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                raise ValueError(f"Invalid GENCODE GTF record at {path}:{number}")
            chrom, _, feature, start_text, end_text, _, strand, _, attributes_text = fields
            audit[f"feature_{feature}"] += 1
            if feature != "polyA_site":
                continue
            start, end = int(start_text), int(end_text)
            if start < 1 or end < start or strand not in {"+", "-"}:
                raise ValueError(f"Invalid GENCODE coordinate at {path}:{number}")
            position = end - 1 if strand == "+" else start - 1
            attributes = _gtf_attributes(attributes_text)
            identifier = attributes.get("transcript_id") or attributes.get("gene_id") or f"GENCODE_{number}"
            gene = attributes.get("gene_id", "")
            key = (chrom, position, strand)
            records.setdefault(key, PasSourceRecord(chrom, position, strand, identifier, gene))
    audit["retained_unique_polyA_site"] = len(records)
    return [records[key] for key in sorted(records)], audit


def _polyadb_member(archive: zipfile.ZipFile, collection: str) -> str:
    suffix = f".PAS.{collection}.tsv"
    matches = [name for name in archive.namelist() if name.endswith(suffix)]
    if len(matches) != 1:
        raise ValueError(f"PolyA_DB archive must contain exactly one *{suffix}; found {matches}")
    return matches[0]


def _parse_polyadb_id(value: str) -> tuple[str, int, str]:
    try:
        chrom, strand, position_text = value.rsplit(":", 2)
        position = int(position_text) - 1
    except (ValueError, TypeError) as exc:
        raise ValueError(f"Invalid PolyA_DB PAS_ID: {value!r}") from exc
    if not chrom or strand not in {"+", "-"} or position < 0:
        raise ValueError(f"Invalid PolyA_DB PAS_ID: {value!r}")
    return chrom, position, strand


def read_polyadb_collection(path: Path, collection: str) -> tuple[list[PasSourceRecord], Counter[str]]:
    if collection not in {"main", "max"}:
        raise ValueError(f"Unsupported PolyA_DB collection: {collection}")
    records: dict[tuple[str, int, str], PasSourceRecord] = {}
    audit: Counter[str] = Counter()
    with zipfile.ZipFile(path) as archive:
        member = _polyadb_member(archive, collection)
        with archive.open(member) as binary:
            import io
            reader = csv.DictReader(io.TextIOWrapper(binary, encoding="utf-8"), delimiter="\t")
            if reader.fieldnames is None or "PAS_ID" not in reader.fieldnames:
                raise ValueError(f"PolyA_DB {member} has no PAS_ID column")
            for row in reader:
                audit["input_rows"] += 1
                chrom, position, strand = _parse_polyadb_id(row["PAS_ID"])
                key = (chrom, position, strand)
                records.setdefault(key, PasSourceRecord(
                    chrom, position, strand, row["PAS_ID"], row.get("GeneSymbol", "")
                ))
    audit["retained_unique_sites"] = len(records)
    audit["duplicate_gene_assignment_rows"] = audit["input_rows"] - len(records)
    return [records[key] for key in sorted(records)], audit


def lift_records(
    records: list[PasSourceRecord], chain: Path, lift_over: str
) -> tuple[list[PasSourceRecord], Counter[str]]:
    audit: Counter[str] = Counter()
    with tempfile.TemporaryDirectory(prefix="rna_ends_polyadb_liftover_") as temporary:
        root = Path(temporary)
        source, mapped, unmapped = root / "source.bed", root / "mapped.bed", root / "unmapped.bed"
        with source.open("w", encoding="utf-8") as handle:
            for index, record in enumerate(records):
                handle.write(
                    f"{record.chrom}\t{record.position}\t{record.position + 1}\tR{index}\t0\t{record.strand}\n"
                )
        subprocess.run(
            [lift_over, "-multiple", str(source), str(chain), str(mapped), str(unmapped)],
            check=True,
        )
        hits: dict[int, list[tuple[str, int, str]]] = {}
        if mapped.is_file():
            with mapped.open(encoding="utf-8") as handle:
                for line in handle:
                    fields = line.rstrip("\n").split("\t")
                    if len(fields) < 6:
                        continue
                    chrom, start, _, token, _, strand = fields[:6]
                    hits.setdefault(int(token[1:]), []).append((chrom, int(start), strand))
        result: list[PasSourceRecord] = []
        for index, record in enumerate(records):
            mapped_hits = hits.get(index, [])
            if len(mapped_hits) != 1:
                audit["unmapped" if not mapped_hits else "multiple"] += 1
                continue
            chrom, position, strand = mapped_hits[0]
            if strand != record.strand:
                audit["strand_changed"] += 1
                continue
            result.append(PasSourceRecord(chrom, position, strand, record.source_id, record.source_gene))
            audit["unique_strand_preserved"] += 1
    return result, audit


def _deduplicate_target_records(
    records: list[PasSourceRecord], audit: Counter[str]
) -> list[PasSourceRecord]:
    unique: dict[tuple[str, int, str], PasSourceRecord] = {}
    for record in records:
        unique.setdefault((record.chrom, record.position, record.strand), record)
    audit["target_coordinate_collisions"] = len(records) - len(unique)
    return [unique[key] for key in sorted(unique)]


def _write_bed(path: Path, records: list[PasSourceRecord]) -> None:
    with gzip.GzipFile(filename=str(path), mode="wb", mtime=0) as compressed:
        import io
        with io.TextIOWrapper(compressed, encoding="utf-8", newline="") as handle:
            for record in sorted(records, key=lambda row: (row.chrom, row.position, row.strand, row.source_id)):
                handle.write(
                    f"{record.chrom}\t{record.position}\t{record.position + 1}\t"
                    f"{record.source_id}\t0\t{record.strand}\t{record.source_gene}\tfalse\n"
                )


def prepare_sources(
    *,
    species: str,
    assembly: str,
    annotation_release: str,
    gencode_polya_gtf: Path,
    polyadb_zip: Path,
    polyadb_assembly: str,
    output: Path,
    download_date: str,
    liftover_chain: Path | None = None,
    lift_over: str = "liftOver",
) -> dict[str, object]:
    profile = SOURCE_PROFILES.get(species)
    if profile is None or (assembly, annotation_release) != (profile["assembly"], profile["annotation"]):
        raise ValueError(f"Unsupported target profile: {species}/{assembly}/{annotation_release}")
    if polyadb_assembly != profile["polyadb_assembly"]:
        raise ValueError(f"Expected PolyA_DB assembly {profile['polyadb_assembly']} for {species}")
    for path in (gencode_polya_gtf, polyadb_zip):
        if not path.is_file():
            raise ValueError(f"Source input does not exist: {path}")
    foreign_polyadb = polyadb_assembly not in {assembly, "hg38" if assembly == "GRCh38" else assembly}
    if foreign_polyadb and (liftover_chain is None or not liftover_chain.is_file()):
        raise ValueError(f"{polyadb_assembly} PolyA_DB requires a lift-over chain to {assembly}")
    output.mkdir(parents=True, exist_ok=True)
    names = {
        "gencode": f"gencode_polyA.{assembly}.bed.gz",
        "main": f"polyadb_v4.1_main.{assembly}.bed.gz",
        "max": f"polyadb_v4.1_max.{assembly}.bed.gz",
    }
    targets = {key: output / name for key, name in names.items()}
    occupied = [str(path) for path in (*targets.values(), output / "source_manifest.json") if path.exists()]
    if occupied:
        raise ValueError(f"Refusing to overwrite prepared PAS sources: {occupied}")

    gencode, gencode_audit = read_gencode_polya_sites(gencode_polya_gtf)
    main, main_audit = read_polyadb_collection(polyadb_zip, "main")
    maximum, max_audit = read_polyadb_collection(polyadb_zip, "max")
    lift_audit: dict[str, dict[str, int]] = {}
    if foreign_polyadb:
        assert liftover_chain is not None
        main, main_lift = lift_records(main, liftover_chain, lift_over)
        maximum, max_lift = lift_records(maximum, liftover_chain, lift_over)
        main = _deduplicate_target_records(main, main_lift)
        maximum = _deduplicate_target_records(maximum, max_lift)
        lift_audit = {"main": dict(main_lift), "max": dict(max_lift)}
    _write_bed(targets["gencode"], gencode)
    _write_bed(targets["main"], main)
    _write_bed(targets["max"], maximum)

    raw_inputs = [gencode_polya_gtf, polyadb_zip] + ([liftover_chain] if liftover_chain else [])
    manifest = {
        "sources": [
            {
                "file": names["gencode"], "url": profile["gencode_url"],
                "release": annotation_release, "download_date": download_date,
                "license": "GENCODE data-use policy; citation requested",
                "sha256": sha256(targets["gencode"]),
                "derived_from": {"file": gencode_polya_gtf.name, "sha256": sha256(gencode_polya_gtf)},
            },
            {
                "file": names["main"], "url": profile["polyadb_url"],
                "release": "PolyA_DB_v4.1_Main", "download_date": download_date,
                "license": "Freely available from PolyA_DB; cite Yu et al. 2025",
                "sha256": sha256(targets["main"]),
                "derived_from": {"file": polyadb_zip.name, "sha256": sha256(polyadb_zip)},
            },
            {
                "file": names["max"], "url": profile["polyadb_url"],
                "release": "PolyA_DB_v4.1_Max", "download_date": download_date,
                "license": "Freely available from PolyA_DB; cite Yu et al. 2025",
                "sha256": sha256(targets["max"]),
                "derived_from": {"file": polyadb_zip.name, "sha256": sha256(polyadb_zip)},
            },
        ]
    }
    (output / "source_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    provenance = {
        "species": species, "assembly": assembly, "annotation_release": annotation_release,
        "polyadb_source_assembly": polyadb_assembly, "coordinate_system": "0-based half-open",
        "gencode_feature_policy": "polyA_site_only",
        "polyasite_policy": "deferred_from_v1",
        "prepared_at": datetime.now(timezone.utc).isoformat(),
        "raw_inputs": [{"path": str(path.resolve()), "sha256": sha256(path)} for path in raw_inputs],
        "coordinate_conversion": None if not foreign_polyadb else {
            "source": polyadb_assembly,
            "target": assembly,
            "chain_url": (
                "https://hgdownload.soe.ucsc.edu/goldenPath/mm10/liftOver/"
                "mm10ToMm39.over.chain.gz"
            ),
            "chain_license": "UCSC Genome Browser EULA",
            "policy": "unique_mapping_and_strand_preserved",
        },
        "audits": {"gencode": dict(gencode_audit), "polyadb_main": dict(main_audit),
                   "polyadb_max": dict(max_audit), "liftover": lift_audit},
        "outputs": {key: {"path": str(path.resolve()), "sha256": sha256(path)}
                    for key, path in targets.items()},
    }
    (output / "preparation_provenance.json").write_text(
        json.dumps(provenance, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    with (output / "SHA256SUMS").open("w", encoding="utf-8") as handle:
        for filename in sorted((*names.values(), "source_manifest.json", "preparation_provenance.json")):
            handle.write(f"{sha256(output / filename)}  {filename}\n")
    return provenance
