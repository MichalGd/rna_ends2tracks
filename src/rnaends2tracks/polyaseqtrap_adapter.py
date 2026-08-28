from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass, field
from pathlib import Path

import pysam

from .apa_a import annotate_site, load_genes, reverse_complement, transcript_end
from .paths import workflow_asset


POLYASEQTRAP_COMMIT = "176ea2884ff1c6be7c64bc44fa7661d82d90e718"
DEEPIP_COMMIT = "988564875d002b6d5d48d8dfb228cba3492dd776"
DEEPIP_MODEL_SHA256 = {
    "human": "d74138c788102ae57a50664b6858a0b79951430fee9bdcc93b07f9b1ba16edf1",
    "mouse": "aba432c85ef6c14e56a6222106acaffbcc3b9131a86508afdf66311fe57123e9",
}


@dataclass
class Candidate:
    chrom: str
    strand: str
    position: int
    sample_counts: dict[str, int] = field(default_factory=dict)
    levels: set[str] = field(default_factory=set)
    repeat_detected: bool = False
    source_sites: int = 0
    candidate_id: str = ""
    deepip_status: str = "not_evaluated"
    deepip_prediction: str = "NA"
    deepip_score: str = "NA"

    @property
    def total(self) -> int:
        return sum(self.sample_counts.values())

    @property
    def samples(self) -> int:
        return sum(value > 0 for value in self.sample_counts.values())


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _load_json(path: Path, label: str) -> dict[str, object]:
    if not path.is_file():
        raise RuntimeError(f"{label} is unavailable: {path}")
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError) as exc:
        raise RuntimeError(f"{label} is invalid JSON: {path}") from exc
    if not isinstance(value, dict):
        raise RuntimeError(f"{label} must contain a JSON object: {path}")
    return value


def verify_manifests(installation: dict[str, object], accepted: dict[str, object], species: str,
                     assembly: str) -> tuple[Path, Path, str]:
    if accepted.get("schema_version") != 1 or accepted.get("status") != "accepted":
        raise RuntimeError("APA-B validation manifest has not been accepted")
    if assembly not in accepted.get("assemblies", []):
        raise RuntimeError(f"APA-B validation manifest does not cover {assembly}")
    if accepted.get("umi_present") is not False or accepted.get("coordinate_deduplication") is not False:
        raise RuntimeError("APA-B validation must record no UMI and no coordinate deduplication")
    if "quantseq_rev_v2_se" not in accepted.get("library_protocols", []):
        raise RuntimeError("APA-B validation does not cover QuantSeq REV V2 single-end libraries")
    pilot = accepted.get("pilot", {})
    if (not isinstance(pilot, dict) or pilot.get("synthetic_pass") is not True
            or pilot.get("real_quantseq_rev_canaries", {}).get(assembly) != "PASS"):
        raise RuntimeError(f"APA-B validation lacks passing synthetic and {assembly} real canaries")
    if accepted.get("quantseq_rev_adaptation") != "genomewide_no_tail_weighted_PAC":
        raise RuntimeError("APA-B validation covers a different QuantSeq REV adaptation")
    if not str(accepted.get("reviewed_by", "")).strip() or not str(accepted.get("accepted_at", "")).strip():
        raise RuntimeError("APA-B validation lacks reviewer or acceptance timestamp")
    model_path, deepip_script, environment_hash = verify_installation(installation, species)
    engine = installation.get("engine", {})
    model = installation.get("models", {}).get(species, {}) if isinstance(installation.get("models"), dict) else {}
    environment = installation.get("environment", {})
    expected_model = accepted.get("models", {}).get(species, accepted.get("model", {})) \
        if isinstance(accepted.get("models", {}), dict) else accepted.get("model", {})
    checks = (
        (engine.get("source_commit"), accepted.get("engine", {}).get("source_commit"), "engine commit"),
        (model.get("sha256"), expected_model.get("sha256"), f"{species} DeepIP model"),
        (environment.get("sha256"), accepted.get("environment", {}).get("sha256"), "environment lock"),
    )
    for observed, expected, label in checks:
        if not expected or str(observed).lower() != str(expected).lower():
            raise RuntimeError(f"Installed APA-B {label} does not match the accepted pilot")
    return model_path, deepip_script, environment_hash


def verify_installation(installation: dict[str, object], species: str) -> tuple[Path, Path, str]:
    """Verify immutable engine/model assets without asserting scientific acceptance."""
    engine = installation.get("engine", {})
    models = installation.get("models", {})
    model = models.get(species, {}) if isinstance(models, dict) else {}
    environment = installation.get("environment", {})
    if engine.get("source_commit") != POLYASEQTRAP_COMMIT:
        raise RuntimeError("Installed PolyAseqTrap commit differs from the workflow pin")
    if installation.get("deepip", {}).get("source_commit") != DEEPIP_COMMIT:
        raise RuntimeError("Installed DeepIP commit differs from the workflow pin")
    model_path = Path(str(model.get("path", "")))
    deepip_script = Path(str(installation.get("deepip", {}).get("script", "")))
    if not model_path.is_file() or sha256(model_path) != DEEPIP_MODEL_SHA256[species]:
        raise RuntimeError(f"Pinned {species} DeepIP model is missing or has the wrong checksum")
    if not deepip_script.is_file():
        raise RuntimeError(f"Pinned DeepIP script is unavailable: {deepip_script}")
    environment_hash = str(environment.get("sha256", ""))
    if not re.fullmatch(r"[0-9a-fA-F]{64}", environment_hash):
        raise RuntimeError("Installed APA-B environment lock checksum is missing or invalid")
    return model_path, deepip_script, environment_hash


def extract_end_counts(source: Path, target: Path) -> dict[str, int]:
    """Collapse C0 alignments to transcript-oriented cleavage endpoints without deduplication."""
    counts: dict[tuple[str, str, int], int] = defaultdict(int)
    records = eligible = duplicate_flagged = soft_clipped = 0
    with pysam.AlignmentFile(source, "rb") as incoming:
        for read in incoming.fetch(until_eof=True):
            records += 1
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            eligible += 1
            duplicate_flagged += int(read.is_duplicate)
            position, strand, clipped = transcript_end(read)
            soft_clipped += int(clipped)
            counts[(read.reference_name, strand, position)] += 1
    target.parent.mkdir(parents=True, exist_ok=True)
    with target.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["chrom", "position", "strand", "count"])
        for (chrom, strand, position), count in sorted(counts.items()):
            writer.writerow([chrom, position, strand, count])
    assigned = sum(counts.values())
    if assigned != eligible:
        raise RuntimeError(f"APA-B endpoint extraction changed eligible record count for {source}: {eligible} != {assigned}")
    return {"records_read": records, "eligible_records": eligible, "records_written": assigned,
            "duplicate_flagged_records_retained": duplicate_flagged,
            "end_soft_clipped_records_included": soft_clipped,
            "distinct_exact_end_coordinates": len(counts)}


def _run(command: list[str], log: Path) -> None:
    log.parent.mkdir(parents=True, exist_ok=True)
    with log.open("w", encoding="utf-8") as handle:
        handle.write("COMMAND: " + subprocess.list2cmdline(command) + "\n")
        handle.flush()
        subprocess.run(command, stdout=handle, stderr=subprocess.STDOUT, check=True)


def _truth(value: object) -> bool:
    return str(value).strip().lower() in {"1", "true", "t", "yes", "y"}


def read_candidates(path: Path, sample_id: str) -> list[Candidate]:
    candidates: list[Candidate] = []
    with path.open(encoding="utf-8-sig", newline="") as handle:
        for line, row in enumerate(csv.DictReader(handle, delimiter="\t"), start=2):
            try:
                count = int(float(row["count"]))
                position = int(row["position"])
            except (KeyError, ValueError) as exc:
                raise RuntimeError(f"Invalid PolyAseqTrap candidate at {path}:{line}") from exc
            if count < 1 or row.get("strand") not in {"+", "-"}:
                continue
            candidates.append(Candidate(
                chrom=row["chrom"], strand=row["strand"], position=position,
                sample_counts={sample_id: count}, levels={row.get("coordinate_level", "unknown")},
                repeat_detected=_truth(row.get("repeat_detected", False)), source_sites=1,
            ))
    return candidates


def cluster_candidates(candidates: list[Candidate], gap: int) -> list[Candidate]:
    result: list[Candidate] = []
    candidates.sort(key=lambda row: (row.chrom, row.strand, row.position))
    current: list[Candidate] = []
    key: tuple[str, str] | None = None
    last = -10**18
    for candidate in candidates:
        candidate_key = (candidate.chrom, candidate.strand)
        if current and (candidate_key != key or candidate.position - last > gap):
            result.append(_merge_cluster(current))
            current = []
        current.append(candidate); key = candidate_key; last = candidate.position
    if current:
        result.append(_merge_cluster(current))
    return result


def _merge_cluster(rows: list[Candidate]) -> Candidate:
    support: dict[int, int] = defaultdict(int)
    sample_counts: dict[str, int] = defaultdict(int)
    levels: set[str] = set()
    for row in rows:
        support[row.position] += row.total
        for sample, count in row.sample_counts.items():
            sample_counts[sample] += count
        levels.update(row.levels)
    # Equal maxima resolve to the transcript-distal coordinate, then genomic coordinate.
    maximum = max(support.values())
    tied = [position for position, value in support.items() if value == maximum]
    position = max(tied) if rows[0].strand == "+" else min(tied)
    return Candidate(rows[0].chrom, rows[0].strand, position, dict(sample_counts), levels,
                     any(row.repeat_detected for row in rows), sum(row.source_sites for row in rows))


def _sequence(fasta: pysam.FastaFile, candidate: Candidate, flank: int = 100) -> str | None:
    try:
        length = fasta.get_reference_length(candidate.chrom)
    except KeyError:
        return None
    start, end = candidate.position - flank + 1, candidate.position + flank + 1
    if start < 0 or end > length:
        return None
    sequence = fasta.fetch(candidate.chrom, start, end).upper()
    return sequence if candidate.strand == "+" else reverse_complement(sequence)


def _a_rich(sequence: str | None) -> bool:
    if not sequence:
        return True
    downstream = sequence[100:120]
    return downstream.count("A") >= 12 or bool(re.search(r"A{6,}", downstream))


def deepip_filter(candidates: list[Candidate], fasta_path: Path, work: Path, script: Path,
                  model: Path) -> tuple[list[Candidate], list[dict[str, object]]]:
    fasta_file = work / "deepip_candidates.fa"
    pending: dict[str, Candidate] = {}
    audit: list[dict[str, object]] = []
    with pysam.FastaFile(fasta_path) as fasta, fasta_file.open("w", encoding="ascii") as handle:
        for index, candidate in enumerate(candidates, start=1):
            identifier = f"candidate_{index:09d}"
            candidate.candidate_id = identifier
            sequence = _sequence(fasta, candidate)
            requires_model = candidate.repeat_detected or _a_rich(sequence)
            if not requires_model:
                candidate.deepip_status = "not_A_rich"
                audit.append({"candidate_id": identifier, "chrom": candidate.chrom, "start": candidate.position,
                              "strand": candidate.strand, "status": "not_A_rich", "prediction": "NA", "retained": True})
            elif sequence is None:
                candidate.deepip_status = "unscorable_reference_edge"
                audit.append({"candidate_id": identifier, "chrom": candidate.chrom, "start": candidate.position,
                              "strand": candidate.strand, "status": "unscorable_reference_edge", "prediction": "NA", "retained": False})
            else:
                pending[identifier] = candidate
                handle.write(f">{identifier}/{(len(pending) - 1) % 2}\n{sequence}\n")
        if len(pending) == 1:
            # DeepIP prints classification metrics and otherwise fails on a one-class
            # dummy label set. This duplicate affects metrics only, never predictions.
            identifier, candidate = next(iter(pending.items()))
            sequence = _sequence(fasta, candidate)
            handle.write(f">deepip_metric_sentinel/{(int(identifier.rsplit('_', 1)[1]) + 1) % 2}\n{sequence}\n")
    predictions: dict[str, tuple[int, str]] = {}
    if pending:
        output = work / "deepip_predictions.csv"
        _run([sys.executable, str(script), "-testSeq", str(fasta_file),
              "-trainedModel", str(model), "-outputFile", str(output)],
             work / "deepip.log")
        predictions = parse_deepip(output, set(pending))
    retained = {id(candidate): candidate for candidate in candidates}
    for identifier, candidate in pending.items():
        label, score = predictions[identifier]
        keep = label == 1
        candidate.deepip_status = "DeepIP_scored"
        candidate.deepip_prediction = str(label)
        candidate.deepip_score = score
        audit.append({"candidate_id": identifier, "chrom": candidate.chrom, "start": candidate.position,
                      "strand": candidate.strand, "status": "DeepIP_scored", "prediction": label,
                      "score": score, "retained": keep})
        if not keep:
            retained.pop(id(candidate), None)
    for row in audit:
        if row["retained"] is False:
            index = int(str(row["candidate_id"]).split("_")[-1]) - 1
            retained.pop(id(candidates[index]), None)
    return list(retained.values()), audit


def parse_deepip(path: Path, expected: set[str]) -> dict[str, tuple[int, str]]:
    predictions: dict[str, tuple[int, str]] = {}
    with path.open(encoding="utf-8-sig", newline="") as handle:
        sample = handle.read(4096); handle.seek(0)
        delimiter = "\t" if sample.count("\t") > sample.count(",") else ","
        reader = csv.DictReader(handle, delimiter=delimiter)
        if not reader.fieldnames:
            raise RuntimeError("DeepIP output has no header")
        for row in reader:
            identifier = next((row[name] for name in reader.fieldnames
                               if name.lower() in {"id", "name", "title", "seqid", "sequence_id"}), "").lstrip(">")
            identifier = re.sub(r"[/_:][01]$", "", identifier)
            label_text = next((row[name] for name in reader.fieldnames if "predict_label" == name.lower()), "")
            score = next((row[name] for name in reader.fieldnames if "score" in name.lower() or "prob" in name.lower()), "NA")
            if identifier in expected and re.fullmatch(r"[01](?:[.]0+)?", str(label_text)):
                predictions[identifier] = (int(float(label_text)), str(score))
    missing = expected.difference(predictions)
    if missing:
        raise RuntimeError("DeepIP output lacks predictions for: " + ", ".join(sorted(missing)[:10]))
    return predictions


def write_outputs(candidates: list[Candidate], samples: list[str], gtf: Path, outdir: Path,
                  deepip_audit: list[dict[str, object]], provenance: dict[str, object]) -> None:
    genes, bins = load_genes(str(gtf))
    outdir.mkdir(parents=True, exist_ok=True)
    candidates.sort(key=lambda row: (row.chrom, row.position, row.strand))
    catalog_fields = ["pas_id", "gene_id", "chrom", "start", "end", "strand", "feature_class",
                      "total_count", "supporting_samples", "polyaseqtrap_coordinate_levels",
                      "source_site_count", "polyaseqtrap_repeat_detected", "deepip_status",
                      "deepip_prediction", "deepip_score", "quantseq_rev_adaptation"]
    with (outdir / "pas_catalog.tsv").open("w", encoding="utf-8", newline="") as catalog_handle, \
            (outdir / "pas_counts.tsv").open("w", encoding="utf-8", newline="") as count_handle:
        catalog = csv.DictWriter(catalog_handle, fieldnames=catalog_fields, delimiter="\t", lineterminator="\n")
        counts = csv.DictWriter(count_handle, fieldnames=["pas_id", *samples], delimiter="\t", lineterminator="\n")
        catalog.writeheader(); counts.writeheader()
        for index, candidate in enumerate(candidates, start=1):
            pas_id = f"APA_B_{provenance['assembly']}_{index:09d}"
            gene_id, feature = annotate_site(candidate.chrom, candidate.strand, candidate.position, genes, bins)
            catalog.writerow({"pas_id": pas_id, "gene_id": gene_id, "chrom": candidate.chrom,
                              "start": candidate.position, "end": candidate.position + 1, "strand": candidate.strand,
                              "feature_class": feature, "total_count": candidate.total,
                              "supporting_samples": candidate.samples,
                              "polyaseqtrap_coordinate_levels": ",".join(sorted(candidate.levels)),
                              "source_site_count": candidate.source_sites,
                              "polyaseqtrap_repeat_detected": str(candidate.repeat_detected).lower(),
                              "deepip_status": candidate.deepip_status,
                              "deepip_prediction": candidate.deepip_prediction,
                              "deepip_score": candidate.deepip_score,
                              "quantseq_rev_adaptation": "genomewide_no_tail_weighted_PAC"})
            counts.writerow({"pas_id": pas_id, **{sample: candidate.sample_counts.get(sample, 0) for sample in samples}})
    fields = ["candidate_id", "chrom", "start", "strand", "status", "prediction", "score", "retained"]
    with (outdir / "deepip_audit.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader(); writer.writerows(deepip_audit)
    (outdir / "engine_provenance.json").write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def execute(args: argparse.Namespace) -> int:
    outdir = Path(args.outdir).resolve()
    work = outdir / "adapter_work"
    work.mkdir(parents=True, exist_ok=True)
    installation = _load_json(Path(args.installation_manifest), "APA-B installation manifest")
    if args.pilot_mode:
        model, deepip_script, environment_hash = verify_installation(installation, args.species)
    else:
        accepted = _load_json(Path(args.validation_manifest), "APA-B validation manifest")
        model, deepip_script, environment_hash = verify_manifests(
            installation, accepted, args.species, args.assembly
        )
    with Path(args.bam_manifest).open(encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    samples = [row["sample_id"] for row in rows]
    if len(samples) != len(set(samples)) or not samples:
        raise RuntimeError("BAM manifest requires unique non-empty sample_id values")
    r_script = workflow_asset("scripts/R/polyaseqtrap_quantseq_rev.R")
    audits: dict[str, dict[str, int]] = {}

    def process(row: dict[str, str]) -> tuple[str, Path]:
        sample = row["sample_id"]
        sample_work = work / sample
        ends = sample_work / "exact_end_counts.tsv"
        audits[sample] = extract_end_counts(Path(row["bam"]), ends)
        output = sample_work / "polyaseqtrap_candidates.tsv"
        _run(["Rscript", str(r_script), "--ends", str(ends), "--fasta", str(Path(args.fasta).resolve()),
              "--output", str(output), "--audit", str(sample_work / "polyaseqtrap_audit.tsv"),
              "--cluster-gap", str(args.cluster_gap)], sample_work / "polyaseqtrap.log")
        return sample, output

    candidate_rows: list[Candidate] = []
    workers = max(1, min(int(args.threads), len(rows)))
    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures = [pool.submit(process, row) for row in rows]
        for future in as_completed(futures):
            sample, path = future.result()
            candidate_rows.extend(read_candidates(path, sample))
    candidates = [row for row in cluster_candidates(candidate_rows, args.cluster_gap)
                  if row.total >= args.min_reads and row.samples >= args.min_samples]
    if not candidates:
        raise RuntimeError("PolyAseqTrap weighted clustering produced no project-supported PAS")
    candidates, deepip_audit = deepip_filter(candidates, Path(args.fasta), work, deepip_script, model)
    if not candidates:
        raise RuntimeError("DeepIP rejected every project-supported APA-B candidate")
    provenance = {
        "adapter": "rna_ends2tracks PolyAseqTrap/DeepIP QuantSeq REV genome-wide adapter v1",
        "assembly": args.assembly, "species": args.species,
        "engine": {"name": "PolyAseqTrap", "source_commit": POLYASEQTRAP_COMMIT},
        "deepip": {"source_commit": DEEPIP_COMMIT},
        "model": {"name": "DeepIP", "sha256": DEEPIP_MODEL_SHA256[args.species]},
        "environment": {"sha256": environment_hash},
        "pilot_mode": bool(args.pilot_mode),
        "umi_present": False, "coordinate_deduplication": False,
        "duplicate_flagged_records_retained": sum(row["duplicate_flagged_records_retained"] for row in audits.values()),
        "input_records": sum(row["records_read"] for row in audits.values()),
        "endpoint_assigned_records": sum(row["records_written"] for row in audits.values()),
        "quantseq_rev_adaptation": (
            "genome-wide transcript 3-prime endpoints; PolyAseqTrap simpleCluster weighted PAC clustering; "
            "species-specific DeepIP; no FindPTA tail-priority call because QuantSeq REV does not retain tails"
        ),
    }
    write_outputs(candidates, samples, Path(args.gtf), outdir, deepip_audit, provenance)
    (outdir / "adapter_audit.json").write_text(json.dumps({"samples": audits, "final_pas": len(candidates)}, indent=2) + "\n", encoding="utf-8")
    if not args.keep_work:
        shutil.rmtree(work)
    return 0


def parser() -> argparse.ArgumentParser:
    value = argparse.ArgumentParser(description="Pinned PolyAseqTrap/DeepIP adapter for QuantSeq REV APA-B")
    value.add_argument("--bam-manifest", required=True)
    value.add_argument("--fasta", required=True)
    value.add_argument("--gtf", required=True)
    value.add_argument("--species", required=True, choices=sorted(DEEPIP_MODEL_SHA256))
    value.add_argument("--assembly", required=True, choices=["GRCh38", "GRCm39"])
    value.add_argument("--outdir", required=True)
    value.add_argument("--threads", type=int, default=1)
    value.add_argument("--validation-manifest", default="")
    value.add_argument(
        "--pilot-mode", action="store_true",
        help="Run a validation canary from a verified installation before an acceptance manifest exists",
    )
    value.add_argument("--installation-manifest", default=os.environ.get("RNA_ENDS2TRACKS_APA_B_INSTALLATION_MANIFEST", ""))
    value.add_argument("--cluster-gap", type=int, default=24)
    value.add_argument("--min-reads", type=int, default=5)
    value.add_argument("--min-samples", type=int, default=2)
    value.add_argument("--keep-work", action="store_true")
    return value


def main() -> None:
    args = parser().parse_args()
    if not args.installation_manifest:
        raise SystemExit("--installation-manifest or RNA_ENDS2TRACKS_APA_B_INSTALLATION_MANIFEST is required")
    if not args.pilot_mode and not args.validation_manifest:
        raise SystemExit("--validation-manifest is required unless --pilot-mode is used")
    raise SystemExit(execute(args))


if __name__ == "__main__":
    main()
