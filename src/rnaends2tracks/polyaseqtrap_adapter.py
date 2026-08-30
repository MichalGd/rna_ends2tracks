from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import time
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
from dataclasses import dataclass, field
from pathlib import Path

import pysam

from .apa_a import annotate_site, load_genes, reverse_complement, transcript_end
from .mcell2019 import is_end_defining_read
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
                     assembly: str, library_protocol: str) -> tuple[Path, Path, str]:
    if accepted.get("schema_version") != 1 or accepted.get("status") != "accepted":
        raise RuntimeError("APA-B validation manifest has not been accepted")
    if assembly not in accepted.get("assemblies", []):
        raise RuntimeError(f"APA-B validation manifest does not cover {assembly}")
    if accepted.get("umi_present") is not False or accepted.get("coordinate_deduplication") is not False:
        raise RuntimeError("APA-B validation must record no UMI and no coordinate deduplication")
    protocols = accepted.get("library_protocols")
    if protocols is None:
        # Schema-v1 manifests accepted before PE support implicitly cover only
        # the original QuantSeq REV V2 single-end protocol.
        protocols = ["quantseq_rev_v2_se"]
    if not isinstance(protocols, list) or library_protocol not in protocols:
        raise RuntimeError(f"APA-B validation does not cover {library_protocol}")
    pilot = accepted.get("pilot", {})
    if (not isinstance(pilot, dict) or pilot.get("synthetic_pass") is not True
            or pilot.get("c1_c1s_reuse_equivalent") is not True
            or pilot.get("real_quantseq_rev_canaries", {}).get(assembly) != "PASS"):
        raise RuntimeError(
            f"APA-B validation lacks endpoint-reuse equivalence, synthetic, or {assembly} real canaries"
        )
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
        (installation.get("workflow_adapter", {}).get("source_commit"),
         accepted.get("workflow_adapter", {}).get("source_commit"), "workflow adapter commit"),
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
    workflow_adapter = installation.get("workflow_adapter", {})
    engine = installation.get("engine", {})
    models = installation.get("models", {})
    model = models.get(species, {}) if isinstance(models, dict) else {}
    environment = installation.get("environment", {})
    if (not isinstance(workflow_adapter, dict)
            or not str(workflow_adapter.get("release", "")).strip()
            or not re.fullmatch(r"[0-9a-fA-F]{7,40}", str(workflow_adapter.get("source_commit", "")))):
        raise RuntimeError("Installed APA-B workflow-adapter release/commit pin is missing or invalid")
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


def extract_end_counts(source: Path, target: Path, library_layout: str = "SE") -> dict[str, int]:
    """Collapse C0 alignments to transcript-oriented cleavage endpoints without deduplication."""
    counts: dict[tuple[str, str, int], int] = defaultdict(int)
    records = eligible = duplicate_flagged = soft_clipped = non_end_defining_mates = 0
    with pysam.AlignmentFile(source, "rb") as incoming:
        for read in incoming.fetch(until_eof=True):
            records += 1
            if not is_end_defining_read(read, library_layout):
                non_end_defining_mates += 1
                continue
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
            "non_end_defining_mate_records": non_end_defining_mates,
            "distinct_exact_end_coordinates": len(counts)}


def _run(command: list[str], log: Path, environment: dict[str, str] | None = None) -> None:
    log.parent.mkdir(parents=True, exist_ok=True)
    with log.open("w", encoding="utf-8") as handle:
        handle.write("COMMAND: " + subprocess.list2cmdline(command) + "\n")
        handle.flush()
        subprocess.run(command, env=environment, stdout=handle, stderr=subprocess.STDOUT, check=True)


def environment_executable(name: str) -> str:
    """Resolve a tool from the same immutable environment as this Python."""
    environment_bin = Path(sys.executable).resolve().parent
    candidates = (environment_bin / name, environment_bin / f"{name}.exe")
    executable = next((candidate for candidate in candidates if candidate.is_file()), None)
    if executable is None:
        raise RuntimeError(
            f"Required APA-B executable is absent beside {sys.executable}: {name}"
        )
    return str(executable)


def r_subprocess_environment() -> dict[str, str]:
    """Remove parent-shell R overrides and prioritize the APA-B environment."""
    environment = os.environ.copy()
    for variable in ("R_HOME", "R_LIBS", "R_LIBS_USER"):
        environment.pop(variable, None)
    environment_bin = str(Path(sys.executable).resolve().parent)
    environment["PATH"] = environment_bin + os.pathsep + environment.get("PATH", "")
    # Parallelism is across samples. Keep each R worker single-threaded so its
    # BLAS/OpenMP runtime cannot silently oversubscribe the global CPU budget.
    environment.update({
        "OMP_NUM_THREADS": "1",
        "MKL_NUM_THREADS": "1",
        "OPENBLAS_NUM_THREADS": "1",
        "BLIS_NUM_THREADS": "1",
    })
    return environment


def _file_signature(paths: list[Path], parameters: dict[str, object]) -> str:
    digest = hashlib.sha256(json.dumps(parameters, sort_keys=True).encode())
    for path in sorted(paths, key=lambda value: str(value)):
        stat = path.stat()
        digest.update(f"{path.resolve()}:{stat.st_size}:{stat.st_mtime_ns}".encode())
    return digest.hexdigest()


def _checkpoint_valid(path: Path, signature: str, outputs: list[Path]) -> dict[str, object] | None:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None
    if payload.get("schema_version") != 1 or payload.get("signature") != signature:
        return None
    records = payload.get("outputs", [])
    if not isinstance(records, list) or len(records) != len(outputs):
        return None
    observed = {str(output.resolve()): output for output in outputs}
    for record in records:
        if not isinstance(record, dict) or str(record.get("path", "")) not in observed:
            return None
        output = observed[str(record["path"])]
        if not output.is_file():
            return None
        stat = output.stat()
        if stat.st_size != record.get("size") or stat.st_mtime_ns != record.get("mtime_ns"):
            return None
    metadata = payload.get("metadata", {})
    return metadata if isinstance(metadata, dict) else {}


def _write_checkpoint(path: Path, signature: str, outputs: list[Path], metadata: dict[str, object]) -> None:
    payload = {
        "schema_version": 1,
        "signature": signature,
        "outputs": [
            {"path": str(output.resolve()), "size": output.stat().st_size,
             "mtime_ns": output.stat().st_mtime_ns}
            for output in outputs
        ],
        "metadata": metadata,
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    temporary.replace(path)


def _validate_source_receipt(receipt_path: Path, required: list[Path]) -> None:
    receipt = _load_json(receipt_path, "exact-end sample receipt")
    if receipt.get("module") != "exact_ends_sample" or receipt.get("exit_status") != 0:
        raise RuntimeError(f"Exact-end receipt is not successful: {receipt_path}")
    records = {
        str(Path(str(record.get("path", ""))).resolve()): record
        for record in receipt.get("outputs", []) if isinstance(record, dict)
    }
    for path in required:
        resolved = str(path.resolve())
        record = records.get(resolved)
        if record is None:
            raise RuntimeError(f"Exact-end receipt does not cover reusable input: {path}")
        stat = path.stat()
        if stat.st_size != record.get("size"):
            raise RuntimeError(f"Reusable exact-end input changed size after receipt: {path}")
        if record.get("validation") == "sha256":
            if sha256(path) != record.get("sha256"):
                raise RuntimeError(f"Reusable exact-end input failed checksum validation: {path}")
        elif stat.st_mtime_ns != record.get("mtime_ns"):
            raise RuntimeError(f"Reusable exact-end input changed after receipt: {path}")


def _read_position_counts(path: Path) -> dict[tuple[str, str, int], int]:
    counts: dict[tuple[str, str, int], int] = defaultdict(int)
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        for line, row in enumerate(csv.DictReader(handle, delimiter="\t"), start=2):
            try:
                chrom = row["chrom"]
                strand = row["strand"]
                position = int(row.get("start", row.get("position", "")))
                count = int(row["count"])
            except (KeyError, ValueError) as exc:
                raise RuntimeError(f"Invalid reusable endpoint at {path}:{line}") from exc
            if strand not in {"+", "-"} or count < 1:
                raise RuntimeError(f"Invalid reusable endpoint at {path}:{line}")
            counts[(chrom, strand, position)] += count
    return counts


def reuse_exact_end_counts(row: dict[str, str], target: Path) -> dict[str, object]:
    """Reconstruct the adapter's raw endpoint universe from validated C1+C1S."""
    c1 = Path(row["c1"])
    c1s = Path(row["c1s"])
    audit_path = Path(row["exact_end_audit"])
    receipt_path = Path(row["exact_end_receipt"])
    required = [c1, c1s, audit_path]
    missing = [str(path) for path in [*required, receipt_path] if not path.is_file()]
    if missing:
        raise RuntimeError("Reusable exact-end inputs are incomplete: " + ", ".join(missing))
    _validate_source_receipt(receipt_path, required)
    audit = _load_json(audit_path, "exact-end sample audit")
    c1_counts = _read_position_counts(c1)
    c1s_counts = _read_position_counts(c1s)
    observed_c1 = sum(c1_counts.values())
    observed_c1s = sum(c1s_counts.values())
    if (int(audit.get("C1", -1)) != observed_c1
            or int(audit.get("C1S", -1)) != observed_c1s
            or int(audit.get("C0", -1)) != observed_c1 + observed_c1s):
        raise RuntimeError(f"C0=C1+C1S reuse invariant failed for {row['sample_id']}")
    combined: dict[tuple[str, str, int], int] = defaultdict(int)
    for source in (c1_counts, c1s_counts):
        for key, count in source.items():
            combined[key] += count
    target.parent.mkdir(parents=True, exist_ok=True)
    with target.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["chrom", "position", "strand", "count"])
        for (chrom, strand, position), count in sorted(combined.items()):
            writer.writerow([chrom, position, strand, count])
    total = observed_c1 + observed_c1s
    return {
        "records_read": total,
        "eligible_records": total,
        "records_written": total,
        "duplicate_flagged_records_retained": int(
            audit.get("duplicate_flagged_C0", audit.get("duplicate_flagged", 0))
        ),
        "duplicate_flagged_count_scope": (
            "C0" if "duplicate_flagged_C0" in audit else "legacy_all_alignment_records"
        ),
        "end_soft_clipped_records_included": observed_c1s,
        "distinct_exact_end_coordinates": len(combined),
        "endpoint_source": "receipt_validated_C1_plus_C1S",
    }


def _prepare_endpoint_job(row: dict[str, str], target: Path, source_mode: str) -> tuple[str, dict[str, object]]:
    reusable_paths = [Path(row.get(key, "")) for key in ("c1", "c1s", "exact_end_audit", "exact_end_receipt")]
    reusable = all(str(path) and path.is_file() for path in reusable_paths)
    partial = any(str(path) and path.is_file() for path in reusable_paths)
    if source_mode == "exact_ends" and not reusable:
        raise RuntimeError(f"APA-B exact-end reuse was required but is unavailable for {row['sample_id']}")
    if source_mode != "bam" and reusable:
        audit = reuse_exact_end_counts(row, target)
    elif source_mode == "auto" and partial:
        raise RuntimeError(f"APA-B found incomplete reusable exact-end inputs for {row['sample_id']}")
    else:
        audit = {
            **extract_end_counts(Path(row["bam"]), target, row.get("library_layout", "SE")),
            "endpoint_source": "bam_fallback",
        }
    return row["sample_id"], audit


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
                  model: Path, threads: int = 1) -> tuple[list[Candidate], list[dict[str, object]]]:
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
        runtime = os.environ.copy()
        runtime.update({
            "OMP_NUM_THREADS": str(max(1, threads)),
            "MKL_NUM_THREADS": str(max(1, threads)),
            "TF_NUM_INTRAOP_THREADS": str(max(1, threads)),
            "TF_NUM_INTEROP_THREADS": str(max(1, min(2, threads))),
        })
        _run([sys.executable, str(script), "-testSeq", str(fasta_file),
              "-trainedModel", str(model), "-outputFile", str(output)],
             work / "deepip.log", environment=runtime)
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
    with Path(args.bam_manifest).open(encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    samples = [row["sample_id"] for row in rows]
    if len(samples) != len(set(samples)) or not samples:
        raise RuntimeError("BAM manifest requires unique non-empty sample_id values")
    layouts = {str(row.get("library_layout", "SE")).upper() for row in rows}
    if layouts not in ({"SE"}, {"PE"}):
        raise RuntimeError("BAM manifest must use one consistent SE or PE library layout")
    expected_layout = "PE" if args.library_protocol.endswith("_pe") else "SE"
    if layouts != {expected_layout}:
        raise RuntimeError(
            f"BAM-manifest layout {sorted(layouts)} does not match {args.library_protocol}"
        )
    if args.pilot_mode:
        model, deepip_script, environment_hash = verify_installation(installation, args.species)
    else:
        accepted = _load_json(Path(args.validation_manifest), "APA-B validation manifest")
        model, deepip_script, environment_hash = verify_manifests(
            installation, accepted, args.species, args.assembly, args.library_protocol
        )
    r_script = workflow_asset("scripts/R/polyaseqtrap_quantseq_rev.R")
    audits: dict[str, dict[str, object]] = {}
    endpoint_paths: dict[str, Path] = {}
    endpoint_jobs: list[tuple[dict[str, str], Path, str, Path]] = []
    endpoint_started = time.monotonic()
    for row in rows:
        sample = row["sample_id"]
        sample_work = work / sample
        ends = sample_work / "exact_end_counts.tsv"
        endpoint_paths[sample] = ends
        if args.endpoint_source == "bam":
            sources = [Path(row["bam"])]
        else:
            reusable = [Path(row.get(key, "")) for key in ("c1", "c1s", "exact_end_audit", "exact_end_receipt")]
            sources = reusable if all(path.is_file() for path in reusable) else [Path(row["bam"])]
        signature = _file_signature(sources, {
            "stage": "endpoint", "sample": sample, "source_mode": args.endpoint_source,
            "adapter": "genomewide_no_tail_weighted_PAC",
        })
        checkpoint = sample_work / "endpoint_checkpoint.json"
        metadata = _checkpoint_valid(checkpoint, signature, [ends])
        if metadata is not None:
            audits[sample] = metadata
            print(f"APA-B endpoint {sample} reused", flush=True)
        else:
            endpoint_jobs.append((row, ends, signature, checkpoint))

    endpoint_workers = max(1, min(int(args.endpoint_workers or args.threads), len(endpoint_jobs) or 1))
    completed_endpoints = len(rows) - len(endpoint_jobs)
    if endpoint_jobs:
        with ProcessPoolExecutor(max_workers=endpoint_workers) as pool:
            futures = {
                pool.submit(_prepare_endpoint_job, row, ends, args.endpoint_source):
                    (row["sample_id"], ends, signature, checkpoint)
                for row, ends, signature, checkpoint in endpoint_jobs
            }
            for future in as_completed(futures):
                sample, ends, signature, checkpoint = futures[future]
                _sample, audit = future.result()
                audits[sample] = audit
                _write_checkpoint(checkpoint, signature, [ends], audit)
                completed_endpoints += 1
                elapsed = time.monotonic() - endpoint_started
                eta = elapsed / completed_endpoints * (len(rows) - completed_endpoints)
                print(
                    f"APA-B endpoint {sample} completed ({completed_endpoints}/{len(rows)}); "
                    f"elapsed={elapsed:.0f}s; ETA~{eta:.0f}s",
                    flush=True,
                )

    def cluster_sample(row: dict[str, str], ends: Path) -> tuple[str, Path, Path, str]:
        sample = row["sample_id"]
        sample_work = work / sample
        output = sample_work / "polyaseqtrap_candidates.tsv"
        audit_output = sample_work / "polyaseqtrap_audit.tsv"
        signature = _file_signature([ends, Path(args.fasta), r_script], {
            "stage": "polyaseqtrap_cluster", "sample": sample, "cluster_gap": args.cluster_gap,
        })
        checkpoint = sample_work / "cluster_checkpoint.json"
        if _checkpoint_valid(checkpoint, signature, [output, audit_output]) is not None:
            return sample, output, checkpoint, signature
        command = [environment_executable("Rscript"), str(r_script), "--ends", str(ends),
                   "--fasta", str(Path(args.fasta).resolve()), "--output", str(output),
                   "--audit", str(audit_output),
                   "--cluster-gap", str(args.cluster_gap)]
        log = sample_work / "polyaseqtrap.log"
        log.parent.mkdir(parents=True, exist_ok=True)
        with log.open("w", encoding="utf-8") as handle:
            handle.write("COMMAND: " + subprocess.list2cmdline(command) + "\n")
            handle.flush()
            subprocess.run(command, env=r_subprocess_environment(), stdout=handle,
                           stderr=subprocess.STDOUT, check=True)
        _write_checkpoint(checkpoint, signature, [output, audit_output], {"sample_id": sample})
        return sample, output, checkpoint, signature

    candidate_rows: list[Candidate] = []
    workers = max(1, min(int(args.cluster_workers or args.threads), len(rows)))
    cluster_started = time.monotonic()
    completed_clusters = 0
    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures = [pool.submit(cluster_sample, row, endpoint_paths[row["sample_id"]]) for row in rows]
        for future in as_completed(futures):
            sample, path, _checkpoint, _signature = future.result()
            candidate_rows.extend(read_candidates(path, sample))
            completed_clusters += 1
            elapsed = time.monotonic() - cluster_started
            eta = elapsed / completed_clusters * (len(rows) - completed_clusters)
            print(
                f"APA-B PolyAseqTrap {sample} completed ({completed_clusters}/{len(rows)}); "
                f"elapsed={elapsed:.0f}s; ETA~{eta:.0f}s",
                flush=True,
            )
    candidates = [row for row in cluster_candidates(candidate_rows, args.cluster_gap)
                  if row.total >= args.min_reads and row.samples >= args.min_samples]
    if not candidates:
        raise RuntimeError("PolyAseqTrap weighted clustering produced no project-supported PAS")
    print(f"APA-B DeepIP started for {len(candidates)} supported candidates", flush=True)
    candidates, deepip_audit = deepip_filter(
        candidates, Path(args.fasta), work, deepip_script, model,
        threads=int(args.deepip_threads or args.threads),
    )
    print(f"APA-B DeepIP completed; {len(candidates)} candidates retained", flush=True)
    if not candidates:
        raise RuntimeError("DeepIP rejected every project-supported APA-B candidate")
    provenance = {
        "adapter": "rna_ends2tracks PolyAseqTrap/DeepIP QuantSeq REV genome-wide adapter v1",
        "workflow_adapter": installation.get("workflow_adapter", {}),
        "assembly": args.assembly, "species": args.species,
        "library_protocol": args.library_protocol,
        "library_layouts": sorted(layouts),
        "engine": {"name": "PolyAseqTrap", "source_commit": POLYASEQTRAP_COMMIT},
        "deepip": {"source_commit": DEEPIP_COMMIT},
        "model": {"name": "DeepIP", "sha256": DEEPIP_MODEL_SHA256[args.species]},
        "environment": {"sha256": environment_hash},
        "pilot_mode": bool(args.pilot_mode),
        "umi_present": False, "coordinate_deduplication": False,
        "duplicate_flagged_records_retained": sum(row["duplicate_flagged_records_retained"] for row in audits.values()),
        "input_records": sum(row["records_read"] for row in audits.values()),
        "endpoint_assigned_records": sum(row["records_written"] for row in audits.values()),
        "endpoint_sources": {sample: row.get("endpoint_source", "unknown") for sample, row in audits.items()},
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
    value.add_argument(
        "--library-protocol", default="quantseq_rev_v2_se",
        choices=["quantseq_rev_v1_se", "quantseq_rev_v2_se", "quantseq_rev_v1_pe", "quantseq_rev_v2_pe"],
    )
    value.add_argument("--outdir", required=True)
    value.add_argument("--threads", type=int, default=1)
    value.add_argument("--endpoint-workers", type=int, default=0,
                       help="Parallel endpoint preparation processes; 0 uses --threads")
    value.add_argument("--cluster-workers", type=int, default=0,
                       help="Parallel PolyAseqTrap sample processes; 0 uses --threads")
    value.add_argument("--deepip-threads", type=int, default=0,
                       help="CPU threads exposed to DeepIP/TensorFlow; 0 uses --threads")
    value.add_argument("--endpoint-source", choices=["auto", "exact_ends", "bam"], default="auto")
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
