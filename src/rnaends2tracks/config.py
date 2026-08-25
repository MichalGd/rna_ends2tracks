from __future__ import annotations

import csv
import gzip
import hashlib
import itertools
import json
import re
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import yaml


class ConfigError(ValueError):
    pass


REQUIRED_COLUMNS = (
    "sample_id",
    "biological_replicate_id",
    "lane_id",
    "fastq_r1",
    "fastq_r2",
    "condition",
    "batch",
    "subject",
    "library_protocol",
    "library_layout",
    "read_length",
    "kit_catalog",
    "umi_present",
)

SUPPORTED_SPECIES = {"human": {"GRCh38"}, "mouse": {"GRCm39"}}
SUPPORTED_PROTOCOLS = {"quantseq_rev_v2_se", "quantseq_rev_v1_se"}
SAFE_ID = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]*$")
SAFE_DESIGN_TERM = re.compile(r"^[A-Za-z][A-Za-z0-9_.]*$")


@dataclass(frozen=True)
class RunPlan:
    project: dict[str, Any]
    samples: list[dict[str, str]]
    sample_rows: list[dict[str, str]]
    contrasts: list[dict[str, Any]]
    reference: dict[str, Any]


def load_yaml(path: str | Path) -> dict[str, Any]:
    path = Path(path).resolve()
    if not path.is_file():
        raise ConfigError(f"Configuration file does not exist: {path}")
    with path.open(encoding="utf-8") as handle:
        value = yaml.safe_load(handle)
    if not isinstance(value, dict):
        raise ConfigError(f"Expected a YAML mapping in {path}")
    value["_config_path"] = str(path)
    return value


def _resolve(value: str, base: Path) -> str:
    if not value:
        return ""
    candidate = Path(value).expanduser()
    if not candidate.is_absolute():
        candidate = base / candidate
    return str(candidate.resolve())


def load_samplesheet(path: str | Path, check_fastqs: bool = True) -> list[dict[str, str]]:
    path = Path(path).resolve()
    if not path.is_file():
        raise ConfigError(f"Samplesheet does not exist: {path}")
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        missing = [column for column in REQUIRED_COLUMNS if column not in (reader.fieldnames or [])]
        if missing:
            raise ConfigError("Samplesheet is missing columns: " + ", ".join(missing))
        rows = [{key: (value or "").strip() for key, value in row.items()} for row in reader]
    if not rows:
        raise ConfigError("Samplesheet contains no data rows")

    base = path.parent
    seen_lanes: set[tuple[str, str]] = set()
    sample_metadata: dict[str, tuple[str, ...]] = {}
    for line_number, row in enumerate(rows, start=2):
        for key in ("sample_id", "biological_replicate_id", "lane_id", "condition", "kit_catalog"):
            if not row[key] or not SAFE_ID.fullmatch(row[key]):
                raise ConfigError(f"Invalid {key!r} on samplesheet line {line_number}: {row[key]!r}")
        lane_key = (row["sample_id"], row["lane_id"])
        if lane_key in seen_lanes:
            raise ConfigError(f"Duplicate sample/lane pair: {lane_key[0]} / {lane_key[1]}")
        seen_lanes.add(lane_key)
        if row["umi_present"].lower() not in {"false", "no", "0"}:
            raise ConfigError(f"UMIs are unsupported; umi_present must be false (line {line_number})")
        layout = row["library_layout"].upper()
        protocol = row["library_protocol"].lower()
        if layout != "SE":
            raise ConfigError(
                f"Only validated QuantSeq REV single-end profiles are enabled; got {layout} on line {line_number}"
            )
        if protocol not in SUPPORTED_PROTOCOLS:
            raise ConfigError(f"Unsupported or unvalidated library_protocol on line {line_number}: {protocol}")
        try:
            if int(row["read_length"]) < 20:
                raise ValueError
        except ValueError as exc:
            raise ConfigError(f"read_length must be an integer >=20 on line {line_number}") from exc
        row["library_layout"] = layout
        row["library_protocol"] = protocol
        row["fastq_r1"] = _resolve(row["fastq_r1"], base)
        row["fastq_r2"] = _resolve(row["fastq_r2"], base)
        if not row["fastq_r1"].endswith((".fastq.gz", ".fq.gz")):
            raise ConfigError(f"FASTQ R1 must be gzip-compressed on line {line_number}: {row['fastq_r1']}")
        if row["fastq_r2"]:
            raise ConfigError(f"fastq_r2 must be empty for the enabled SE profiles (line {line_number})")
        if check_fastqs and not Path(row["fastq_r1"]).is_file():
            raise ConfigError(f"FASTQ R1 does not exist: {row['fastq_r1']}")
        if check_fastqs:
            try:
                with gzip.open(row["fastq_r1"], "rt", encoding="ascii", errors="replace") as fastq:
                    first = fastq.readline(); fastq.readline(); third = fastq.readline()
                if not first.startswith("@") or not third.startswith("+"):
                    raise ConfigError(f"Invalid FASTQ structure: {row['fastq_r1']}")
            except (OSError, EOFError) as exc:
                raise ConfigError(f"Cannot read gzip FASTQ: {row['fastq_r1']}: {exc}") from exc

        invariant_keys = (
            "biological_replicate_id", "condition", "batch", "subject",
            "library_protocol", "library_layout", "read_length", "kit_catalog", "umi_present",
        )
        invariant = tuple(row[key] for key in invariant_keys)
        prior = sample_metadata.setdefault(row["sample_id"], invariant)
        if prior != invariant:
            raise ConfigError(f"Biological metadata differs among lanes for sample {row['sample_id']}")
    return rows


def collapse_samples(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    result: list[dict[str, str]] = []
    seen: set[str] = set()
    for row in rows:
        if row["sample_id"] in seen:
            continue
        seen.add(row["sample_id"])
        result.append({key: value for key, value in row.items() if key not in {"lane_id", "fastq_r1", "fastq_r2"}})
    return result


def generate_contrasts(samples: list[dict[str, str]], order: list[str]) -> list[dict[str, Any]]:
    by_condition: dict[str, set[str]] = {}
    for sample in samples:
        by_condition.setdefault(sample["condition"], set()).add(sample["biological_replicate_id"])
    unknown = sorted(set(by_condition) - set(order))
    absent = sorted(set(order) - set(by_condition))
    if unknown:
        raise ConfigError("condition_order is missing observed conditions: " + ", ".join(unknown))
    if absent:
        raise ConfigError("condition_order contains conditions absent from samples: " + ", ".join(absent))
    eligible = [condition for condition in order if len(by_condition[condition]) >= 2]
    if len(eligible) < 2:
        raise ConfigError("At least two conditions with >=2 biological replicates are required")
    contrasts: list[dict[str, Any]] = []
    for denominator, numerator in itertools.combinations(eligible, 2):
        n_num = len(by_condition[numerator])
        n_den = len(by_condition[denominator])
        pairing_status, n_pairs = _pairing_state(samples, numerator, denominator, "subject")
        paired = pairing_status == "complete"
        contrasts.append({
            "contrast_id": f"{numerator}_vs_{denominator}",
            "factor": "condition",
            "numerator": numerator,
            "denominator": denominator,
            "n_num": n_num,
            "n_den": n_den,
            "paired": paired,
            "n_pairs": n_pairs,
            "pairing_status": pairing_status,
            "design_status": "LOW_REPLICATION_N2" if min(n_num, n_den) == 2 else "valid",
        })
    return contrasts


def _pairing_state(
    samples: list[dict[str, str]], numerator: str, denominator: str, subject_column: str,
) -> tuple[str, int]:
    numerator_values = [sample.get(subject_column, "").strip() for sample in samples if sample["condition"] == numerator]
    denominator_values = [sample.get(subject_column, "").strip() for sample in samples if sample["condition"] == denominator]
    numerator_nonempty = [value for value in numerator_values if value]
    denominator_nonempty = [value for value in denominator_values if value]
    if not numerator_nonempty and not denominator_nonempty:
        return "no_subjects", 0
    if len(numerator_nonempty) != len(numerator_values) or len(denominator_nonempty) != len(denominator_values):
        return "incomplete", 0
    numerator_counts = Counter(numerator_nonempty)
    denominator_counts = Counter(denominator_nonempty)
    if set(numerator_counts).isdisjoint(denominator_counts):
        return "disjoint_subjects", 0
    if numerator_counts == denominator_counts and all(count == 1 for count in numerator_counts.values()):
        return "complete", len(numerator_counts)
    return "incomplete", 0


def resolve_contrast_designs(
    samples: list[dict[str, str]], contrasts: list[dict[str, Any]], project: dict[str, Any],
) -> list[dict[str, Any]]:
    statistics = project.get("statistics", {})
    if not isinstance(statistics, dict):
        raise ConfigError("statistics must be a mapping")
    pairing = statistics.get("pairing", {})
    if not isinstance(pairing, dict):
        raise ConfigError("statistics.pairing must be a mapping")
    mode = str(pairing.get("mode", "none")).lower()
    if mode not in {"none", "auto", "required"}:
        raise ConfigError("statistics.pairing.mode must be one of: none, auto, required")
    subject_column = str(pairing.get("subject_column", "subject")).strip()
    if not SAFE_DESIGN_TERM.fullmatch(subject_column):
        raise ConfigError("statistics.pairing.subject_column must be a safe column name")
    paired_design = str(pairing.get("paired_design", f"~ {subject_column} + condition")).strip()
    incomplete_action = str(pairing.get("incomplete_pair_action", "error")).lower()
    if incomplete_action not in {"error", "unpaired"}:
        raise ConfigError("statistics.pairing.incomplete_pair_action must be 'error' or 'unpaired'")
    default_design = str(project["design"])
    if mode != "none":
        variables = _design_variables(paired_design)
        if "condition" not in variables or subject_column not in variables:
            raise ConfigError(
                "statistics.pairing.paired_design must contain condition and the configured subject column"
            )
        if any(subject_column not in sample for sample in samples):
            raise ConfigError(f"Pairing column is absent from samplesheet: {subject_column}")

    resolved: list[dict[str, Any]] = []
    for original in contrasts:
        contrast = dict(original)
        status, n_pairs = _pairing_state(
            samples, str(contrast["numerator"]), str(contrast["denominator"]), subject_column,
        )
        metadata_paired = status == "complete"
        design_mode = "unpaired"
        formula = default_design
        if mode == "required" and not metadata_paired:
            raise ConfigError(
                f"Complete pairing is required for {contrast['contrast_id']}; observed pairing status: {status}"
            )
        if mode in {"auto", "required"} and metadata_paired:
            design_mode = "paired"
            formula = paired_design
        elif mode == "auto" and status == "incomplete" and incomplete_action == "error":
            raise ConfigError(
                f"Incomplete subject pairing for {contrast['contrast_id']}; fix subject IDs or set "
                "statistics.pairing.incomplete_pair_action: unpaired"
            )
        subset = [
            sample for sample in samples
            if sample["condition"] in {contrast["numerator"], contrast["denominator"]}
        ]
        validate_design(subset, formula)
        contrast.update({
            "paired": metadata_paired,
            "n_pairs": n_pairs if metadata_paired else 0,
            "pairing_status": status,
            "design_mode": design_mode,
            "resolved_design": formula,
            "pairing_column": subject_column if design_mode == "paired" else "",
        })
        resolved.append(contrast)
    return resolved


def _design_variables(design: str) -> list[str]:
    if any(operator in design for operator in ("*", ":", "/", "(")):
        raise ConfigError("Initial implementation supports additive designs only, e.g. '~ batch + condition'")
    right = design.partition("~")[2]
    variables = [item.strip() for item in right.split("+") if item.strip() and item.strip() != "1"]
    if not variables:
        raise ConfigError("Design formula contains no variables")
    for variable in variables:
        if not SAFE_DESIGN_TERM.fullmatch(variable):
            raise ConfigError(f"Unsupported design term: {variable!r}")
    return variables


def _matrix_rank(matrix: list[list[float]], tolerance: float = 1e-10) -> int:
    work = [row[:] for row in matrix]
    rows, columns = len(work), len(work[0]) if work else 0
    rank = 0
    for column in range(columns):
        pivot = next((row for row in range(rank, rows) if abs(work[row][column]) > tolerance), None)
        if pivot is None:
            continue
        work[rank], work[pivot] = work[pivot], work[rank]
        divisor = work[rank][column]
        work[rank] = [value / divisor for value in work[rank]]
        for row in range(rows):
            if row != rank and abs(work[row][column]) > tolerance:
                factor = work[row][column]
                work[row] = [a - factor * b for a, b in zip(work[row], work[rank])]
        rank += 1
    return rank


def validate_design(samples: list[dict[str, str]], design: str) -> None:
    variables = _design_variables(design)
    columns: list[list[float]] = [[1.0] for _ in samples]
    column_count = 1
    for variable in variables:
        if variable not in samples[0]:
            raise ConfigError(f"Design variable is absent from samplesheet: {variable}")
        values = [sample.get(variable, "") for sample in samples]
        if any(value == "" for value in values):
            raise ConfigError(f"Design variable has missing values: {variable}")
        levels = sorted(set(values))
        if len(levels) < 2:
            raise ConfigError(f"Design variable has fewer than two levels: {variable}")
        for level in levels[1:]:
            for row, value in zip(columns, values):
                row.append(float(value == level))
            column_count += 1
    if len(samples) < column_count or _matrix_rank(columns) < column_count:
        raise ConfigError("Design matrix is not full rank; condition may be confounded with batch/subject")


def load_reference(project: dict[str, Any], check_files: bool = True) -> dict[str, Any]:
    ref_value = project.get("reference", {})
    if not isinstance(ref_value, dict):
        raise ConfigError("reference must be a mapping")
    project_base = Path(project["_config_path"]).parent
    if "manifest" in ref_value:
        manifest_path = Path(_resolve(str(ref_value["manifest"]), project_base))
        ref = load_yaml(manifest_path)
        base = manifest_path.parent
    else:
        ref = dict(ref_value)
        base = project_base
    species = str(ref.get("species", "")).lower()
    assembly = str(ref.get("assembly", ""))
    if species not in SUPPORTED_SPECIES:
        raise ConfigError("reference species must be 'human' or 'mouse'")
    if assembly not in SUPPORTED_SPECIES[species]:
        allowed = ", ".join(sorted(SUPPORTED_SPECIES[species]))
        raise ConfigError(f"Unsupported {species} assembly {assembly!r}; validated assemblies: {allowed}")
    for key in ("fasta", "gtf", "star_index", "chrom_sizes"):
        if not ref.get(key):
            raise ConfigError(f"Reference manifest is missing {key}")
        ref[key] = _resolve(str(ref[key]), base)
        if check_files and not Path(ref[key]).exists():
            raise ConfigError(f"Reference asset does not exist ({key}): {ref[key]}")
    if ref.get("pas_atlas"):
        ref["pas_atlas"] = _resolve(str(ref["pas_atlas"]), base)
        if check_files and not Path(ref["pas_atlas"]).is_file():
            raise ConfigError(f"PAS atlas does not exist: {ref['pas_atlas']}")
    if check_files:
        _validate_reference_assets(ref)
    warnings: list[dict[str, str]] = []
    fasta_hash = str(ref.get("fasta_sha256", "")).lower()
    gtf_hash = str(ref.get("gtf_sha256", "")).lower()
    index_fasta_hash = str(ref.get("star_index_fasta_sha256", "")).lower()
    index_gtf_hash = str(ref.get("star_index_gtf_sha256", "")).lower()
    if fasta_hash and index_fasta_hash and fasta_hash != index_fasta_hash:
        raise ConfigError("STAR index FASTA checksum does not match the selected FASTA")
    if gtf_hash and index_gtf_hash and gtf_hash != index_gtf_hash:
        raise ConfigError("STAR index GTF checksum does not match the selected GTF")
    if not (fasta_hash and gtf_hash and index_fasta_hash and index_gtf_hash):
        warnings.append({
            "warning_code": "STAR_INDEX_PROVENANCE_UNVERIFIED",
            "message": "Structural compatibility can be checked, but FASTA/GTF checksums linking the STAR index to the selected references are incomplete.",
        })
    ref["_warnings"] = warnings
    ref.pop("_config_path", None)
    return ref


def _validate_reference_assets(ref: dict[str, Any]) -> None:
    fasta = Path(ref["fasta"])
    if not fasta.is_file() or fasta.suffix == ".gz":
        raise ConfigError("Reference FASTA must be an uncompressed regular file")
    fai = Path(str(fasta) + ".fai")
    if not fai.is_file():
        raise ConfigError(f"Indexed FASTA is required: {fai}")
    fasta_contigs = {line.split("\t", 1)[0] for line in fai.read_text(encoding="utf-8").splitlines() if line}
    with Path(ref["chrom_sizes"]).open(encoding="utf-8") as handle:
        size_contigs = {line.split("\t", 1)[0] for line in handle if line.strip()}
    if fasta_contigs != size_contigs:
        raise ConfigError("FASTA index and chromosome-size contig sets differ")
    star_dir = Path(ref["star_index"])
    for asset in ("Genome", "SA", "SAindex", "chrName.txt", "chrLength.txt", "genomeParameters.txt"):
        if not (star_dir / asset).is_file():
            raise ConfigError(f"STAR index is incomplete; missing {star_dir / asset}")
    star_contigs = set((star_dir / "chrName.txt").read_text(encoding="utf-8").split())
    if star_contigs != fasta_contigs:
        raise ConfigError("STAR index and FASTA contig sets differ")
    star_names = (star_dir / "chrName.txt").read_text(encoding="utf-8").split()
    star_lengths_text = (star_dir / "chrLength.txt").read_text(encoding="utf-8").split()
    if len(star_names) != len(star_lengths_text):
        raise ConfigError("STAR chrName.txt and chrLength.txt have different row counts")
    try:
        star_lengths = {name: int(length) for name, length in zip(star_names, star_lengths_text)}
        fasta_lengths = {
            fields[0]: int(fields[1])
            for line in fai.read_text(encoding="utf-8").splitlines()
            if (fields := line.split("\t")) and len(fields) >= 2
        }
    except ValueError as exc:
        raise ConfigError("Invalid numeric contig length in FASTA or STAR index metadata") from exc
    if star_lengths != fasta_lengths:
        raise ConfigError("STAR index and FASTA contig lengths differ")
    parameters: dict[str, str] = {}
    for line in (star_dir / "genomeParameters.txt").read_text(encoding="utf-8").splitlines():
        key, _, value = line.strip().partition("\t")
        if not value:
            key, _, value = line.strip().partition(" ")
        if key and value:
            parameters[key] = value.strip()
    ref["star_index_version"] = parameters.get("versionGenome", "unknown")
    if parameters.get("sjdbOverhang", "").isdigit():
        ref["star_sjdb_overhang"] = int(parameters["sjdbOverhang"])
    opener = gzip.open if str(ref["gtf"]).endswith(".gz") else open
    gtf_contigs: set[str] = set()
    with opener(ref["gtf"], "rt", encoding="utf-8") as handle:
        for line in handle:
            if line and not line.startswith("#"):
                gtf_contigs.add(line.split("\t", 1)[0])
    unknown = sorted(gtf_contigs - fasta_contigs)
    if unknown:
        raise ConfigError("GTF has contigs absent from FASTA: " + ", ".join(unknown[:10]))
    if ref.get("pas_atlas"):
        opener = gzip.open if str(ref["pas_atlas"]).endswith(".gz") else open
        atlas_contigs: set[str] = set()
        with opener(ref["pas_atlas"], "rt", encoding="utf-8") as handle:
            for line in handle:
                if line and not line.startswith(("#", "track", "browser")):
                    atlas_contigs.add(line.split("\t", 1)[0])
        unknown_atlas = sorted(atlas_contigs - fasta_contigs)
        if unknown_atlas:
            raise ConfigError("PAS atlas has contigs absent from FASTA: " + ", ".join(unknown_atlas[:10]))


def build_plan(config_path: str | Path, samplesheet_path: str | Path, check_inputs: bool = True) -> RunPlan:
    project = load_yaml(config_path)
    if not SAFE_ID.fullmatch(str(project.get("project_id", ""))):
        raise ConfigError("project_id must be a filesystem-safe identifier")
    if project.get("protocol", {}).get("has_umi") not in (False, "false", 0):
        raise ConfigError("protocol.has_umi must be false")
    if project.get("protocol", {}).get("retain_duplicate_flagged_reads") is not True:
        raise ConfigError("protocol.retain_duplicate_flagged_reads must be true for no-UMI QuantSeq")
    profile = str(project.get("protocol", {}).get("profile", "")).lower()
    if profile not in SUPPORTED_PROTOCOLS:
        raise ConfigError(f"Unsupported or unvalidated protocol profile: {profile}")
    design = str(project.get("design", "~ condition"))
    if "condition" not in design:
        raise ConfigError("The model design must include condition")
    rows = load_samplesheet(samplesheet_path, check_fastqs=check_inputs)
    if {row["library_protocol"] for row in rows} != {profile}:
        raise ConfigError("Samplesheet library_protocol does not match project protocol.profile")
    samples = collapse_samples(rows)
    project["design"] = design
    replicate_owner: dict[str, str] = {}
    for sample in samples:
        prior = replicate_owner.setdefault(sample["biological_replicate_id"], sample["sample_id"])
        if prior != sample["sample_id"]:
            raise ConfigError(
                f"biological_replicate_id {sample['biological_replicate_id']} belongs to multiple sample IDs; "
                "technical lanes must share one sample_id"
            )
    validate_design(samples, design)
    order = project.get("condition_order")
    if not isinstance(order, list) or not all(isinstance(x, str) for x in order):
        raise ConfigError("condition_order must be an explicit YAML list")
    if len(order) != len(set(order)):
        raise ConfigError("condition_order must not contain duplicates")
    contrasts = resolve_contrast_designs(samples, generate_contrasts(samples, order), project)
    reference = load_reference(project, check_files=check_inputs)
    expected_overhang = max(int(row["read_length"]) for row in rows) - 1
    observed_overhang = reference.get("star_sjdb_overhang")
    if observed_overhang is not None and observed_overhang != expected_overhang:
        reference["star_overhang_review"] = (
            f"STAR sjdbOverhang={observed_overhang}; declared maximum read length suggests {expected_overhang}. "
            "Reuse is permitted only after reviewing the splice-alignment impact."
        )
        reference.setdefault("_warnings", []).append({
            "warning_code": "STAR_SJDB_OVERHANG_REVIEW",
            "message": reference["star_overhang_review"],
        })
    return RunPlan(project, samples, rows, contrasts, reference)


def write_plan(plan: RunPlan, outdir: str | Path) -> None:
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    with (outdir / "resolved_config.json").open("w", encoding="utf-8") as handle:
        payload = dict(plan.project)
        payload.pop("_config_path", None)
        payload["reference"] = plan.reference
        json.dump(payload, handle, indent=2, sort_keys=True)
        handle.write("\n")
    _write_tsv(outdir / "validated_samples.tsv", plan.samples)
    _write_tsv(outdir / "validated_lanes.tsv", plan.sample_rows)
    _write_tsv(outdir / "contrasts.tsv", plan.contrasts)
    _write_tsv(outdir / "warnings.tsv", plan.reference.get("_warnings", []))


def _write_tsv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def signature_for(paths: list[str | Path], parameters: dict[str, Any]) -> str:
    digest = hashlib.sha256(json.dumps(parameters, sort_keys=True, default=str).encode())
    for raw_path in sorted(map(str, paths)):
        path = Path(raw_path)
        stat = path.stat()
        digest.update(f"{path.resolve()}:{stat.st_size}:{stat.st_mtime_ns}".encode())
    return digest.hexdigest()
