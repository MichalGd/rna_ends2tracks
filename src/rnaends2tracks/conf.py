from __future__ import annotations

import re
import shlex
from pathlib import Path
from typing import Any


class ConfError(ValueError):
    pass


ASSIGNMENT = re.compile(r"^([A-Z][A-Z0-9_]*)[ \t]*=[ \t]*(.*)$")


DEFAULTS: dict[str, str] = {
    "TMP_DIR": "",
    "RUN_GENE_EXPRESSION": "true",
    "RUN_APA_A_MCELL2019": "true",
    # Preserve legacy configurations that predate explicit APA-B settings.
    # New projects copy config/config.conf, where both APA methods are enabled
    # for the audited site scope and guarded by an accepted manifest.
    "RUN_APA_B": "false",
    "APA_B_PILOT_ACCEPTED": "false",
    "APA_B_COMMAND_TEMPLATE": "auto",
    "APA_B_VALIDATION_MANIFEST": "",
    "APA_B_INSTALLATION_MANIFEST": "/opt/conda_envs/rna_ends2tracks-apa-b-0.1.0a10.post6/installation_manifest.json",
    "APA_B_THREADS": "8",
    "APA_B_ENDPOINT_SOURCE": "auto",
    "APA_B_ENDPOINT_PARALLEL_JOBS": "8",
    "APA_B_CLUSTER_PARALLEL_JOBS": "8",
    "APA_B_DEEPIP_THREADS": "8",
    "RUN_DGE_ENRICHMENT": "true",
    "RUN_APA_ENRICHMENT": "true",
    "ENRICHMENT_ORA": "true",
    "ENRICHMENT_GSEA": "true",
    "ENRICHMENT_GO": "true",
    "ENRICHMENT_REACTOME": "true",
    "ENRICHMENT_HALLMARKS": "true",
    "ENRICHMENT_KEGG": "true",
    "ENRICHMENT_RICH_PLOTS": "true",
    "ENRICHMENT_NETWORK_MAX_TERMS": "8",
    "ENRICHMENT_NETWORK_MAX_GENES": "50",
    "ENRICHMENT_PADJ": "0.05",
    "ENRICHMENT_DGE_MIN_ABS_LFC": "1.0",
    "ENRICHMENT_APA_MIN_ABS_DELTA_PAU": "0.10",
    "ENRICHMENT_MIN_GENESET_SIZE": "10",
    "ENRICHMENT_MAX_GENESET_SIZE": "500",
    "ENRICHMENT_PARALLEL_JOBS": "6",
    "RUN_TRACKS": "true",
    "GENERATE_EARLY_C0_TRACKS": "true",
    "RUN_RSEQC": "true",
    "RSEQC_INFER_EXPERIMENT": "true",
    "RSEQC_READ_DISTRIBUTION": "true",
    "RSEQC_GENE_BODY_COVERAGE": "true",
    "RSEQC_MULTIQC": "true",
    "RSEQC_SAMPLE_READS": "200000",
    "RSEQC_MIN_TRANSCRIPT_LENGTH": "100",
    "RSEQC_PARALLEL_JOBS": "6",
    "RSEQC_MEMORY_GB": "4",
    "RUN_FASTQ_SCREEN": "true",
    "FASTQ_SCREEN_CONFIG": "",
    "FASTQ_SCREEN_MISSING_ACTION": "warn",
    "FASTQ_SCREEN_SUBSET": "200000",
    "FASTQ_SCREEN_THREADS": "4",
    "FASTQ_SCREEN_PARALLEL_JOBS": "4",
    "FASTQ_SCREEN_MEMORY_GB": "4",
    "LIBRARY_PROTOCOL": "quantseq_rev_v2_se",
    "LIBRARY_LAYOUT": "single_end",
    "UMI_PRESENT": "false",
    "MAPPING_POLICY": "unique_primary",
    "END_SOFT_CLIP_POLICY": "exclude_and_report",
    "BBDUK_REFERENCE": "adapters,polyA_T",
    "TRIM_QUALITY": "10",
    "MINIMUM_READ_LENGTH": "20",
    "PE_R2_TRIM_5P": "12",
    "ORIENTATION_MIN_FRACTION": "0.75",
    "MIN_REPLICATES_PER_CONDITION": "2",
    "CONDITION_ORDER": "",
    "DESIGN": "~ condition",
    "PAIRING_MODE": "auto",
    "PAIRING_COLUMN": "subject",
    "INCOMPLETE_PAIR_ACTION": "error",
    "AMBIGUOUS_GENE_POLICY": "exclude_statistics",
    "INTERNAL_PRIMING_CONSECUTIVE_BASES": "6",
    "INTERNAL_PRIMING_WINDOW_NT": "10",
    "INTERNAL_PRIMING_MIN_BASES_IN_WINDOW": "7",
    "PAS_MASK_RESCUE_TIER": "core",
    "PAS_DISCOVERY_WINDOW_NT": "30",
    "PAS_DISCOVERY_THRESHOLD": "30",
    "PAS_DISCOVERY_THRESHOLD_OPERATOR": "greater_than",
    "PAS_DISCOVERY_ROUNDS": "2",
    "GENE_DOWNSTREAM_EXTENSION_NT": "6000",
    "FDR": "0.05",
    "MIN_ABS_DELTA_PAU": "0.10",
    "MAX_TOTAL_THREADS": "48",
    "MAX_TOTAL_MEMORY_GB": "384",
    "PARALLEL_DOWNSTREAM_MODULES": "true",
    "DOWNSTREAM_MODULE_PARALLEL_JOBS": "3",
    "PREPROCESS_PARALLEL_JOBS": "4",
    "FASTQC_THREADS": "4",
    "BBDUK_THREADS": "8",
    "BBDUK_MEMORY_GB": "8",
    "STAR_PARALLEL_JOBS": "4",
    "STAR_THREADS": "12",
    "STAR_MEMORY_GB": "48",
    "SAMTOOLS_THREADS": "6",
    "SAMTOOLS_SORT_MEMORY_PER_THREAD_GB": "2",
    "SAMPLE_MERGE_PARALLEL_JOBS": "4",
    "END_EXTRACTION_PARALLEL_JOBS": "8",
    "TRACK_PARALLEL_JOBS": "8",
    "TRACK_THREADS": "4",
    "DGE_CONTRAST_PARALLEL_JOBS": "3",
    "APA_CONTRAST_PARALLEL_JOBS": "4",
    # Empty values preserve the legacy shared APA_CONTRAST_PARALLEL_JOBS
    # setting in existing project configurations.
    "APA_A_CONTRAST_PARALLEL_JOBS": "",
    "APA_B_CONTRAST_PARALLEL_JOBS": "",
    "GENERATE_ALL_READ_TRACKS": "true",
    "GENERATE_EXACT_END_TRACKS": "true",
    "GENERATE_FILTERED_END_TRACKS": "true",
    "GENERATE_REJECTED_END_TRACKS": "true",
    "GENERATE_ACTIVE_PAS_TRACKS": "true",
    "GENERATE_RAW_TRACKS": "true",
    "GENERATE_CPM_TRACKS": "true",
    "GENERATE_DESEQ2_FINAL_TRACKS": "true",
    "GENERATE_DESEQ2_ROBUST_CPM_FINAL_TRACKS": "true",
    "GENERATE_BIGWIGS": "true",
    "RETAIN_BEDGRAPH": "false",
    "UCSC_BIGDATA_URL_PREFIX": "",
    "UCSC_NEGATE_MINUS_TRACKS": "true",
    "UCSC_VIEW_LIMITS": "0:12",
    "CLEANUP_INTERMEDIATES": "true",
    "KEEP_TRIMMED_FASTQ": "false",
    "KEEP_LANE_BAMS": "false",
    "KEEP_APA_SAMPLE_EXTRACTION": "false",
    "KEEP_TRACK_STRAND_BAMS": "false",
    "KEEP_BEDGRAPHS": "false",
}


REQUIRED = {"PROJECT_ID", "SAMPLESHEET", "OUTPUT_DIR"}
REFERENCE_KEYS = {
    "HG38_STAR_INDEX", "HG38_FASTA", "HG38_GTF", "HG38_CHROM_SIZES", "HG38_PAS_ATLAS",
    "HG38_RSEQC_BED",
    "MM39_STAR_INDEX", "MM39_FASTA", "MM39_GTF", "MM39_CHROM_SIZES", "MM39_PAS_ATLAS",
    "MM39_RSEQC_BED",
}
ALLOWED = set(DEFAULTS) | REQUIRED | REFERENCE_KEYS


def _parse_value(raw: str, line_number: int) -> str:
    if any(token in raw for token in ("`", "$(", "${", ";", "&&", "||", "<", ">")):
        raise ConfError(f"Executable shell syntax is forbidden on config line {line_number}")
    try:
        lexer = shlex.shlex(raw, posix=False)
        lexer.whitespace_split = True
        lexer.commenters = "#"
        tokens = list(lexer)
    except ValueError as exc:
        raise ConfError(f"Invalid quoted value on config line {line_number}: {exc}") from exc
    if len(tokens) > 1:
        raise ConfError(f"Config line {line_number} must contain one quoted or unquoted value")
    if not tokens:
        return ""
    value = tokens[0]
    if len(value) >= 2 and value[0] == value[-1] and value[0] in {'"', "'"}:
        value = value[1:-1]
    return value


def read_conf(path: str | Path) -> dict[str, str]:
    source = Path(path).expanduser().resolve()
    if not source.is_file():
        raise ConfError(f"Configuration file does not exist: {source}")
    parsed: dict[str, str] = {}
    for line_number, text in enumerate(source.read_text(encoding="utf-8-sig").splitlines(), start=1):
        stripped = text.strip()
        if not stripped or stripped.startswith("#"):
            continue
        match = ASSIGNMENT.fullmatch(stripped)
        if not match:
            raise ConfError(f"Expected uppercase KEY=value assignment on config line {line_number}")
        key, raw = match.groups()
        if key not in ALLOWED:
            raise ConfError(f"Unknown configuration key on line {line_number}: {key}")
        if key in parsed:
            raise ConfError(f"Duplicate configuration key on line {line_number}: {key}")
        parsed[key] = _parse_value(raw, line_number)
    missing = sorted(key for key in REQUIRED if not parsed.get(key))
    if missing:
        raise ConfError("Missing required configuration keys: " + ", ".join(missing))
    return {**DEFAULTS, **parsed, "_CONFIG_PATH": str(source)}


def _bool(values: dict[str, str], key: str) -> bool:
    value = values[key].strip().lower()
    if value not in {"true", "false"}:
        raise ConfError(f"{key} must be true or false")
    return value == "true"


def _int(values: dict[str, str], key: str, minimum: int = 1) -> int:
    try:
        value = int(values[key])
    except ValueError as exc:
        raise ConfError(f"{key} must be an integer") from exc
    if value < minimum:
        raise ConfError(f"{key} must be >= {minimum}")
    return value


def _int_or(values: dict[str, str], key: str, fallback: str, minimum: int = 1) -> int:
    return _int(values, key if values[key].strip() else fallback, minimum)


def _float(values: dict[str, str], key: str, minimum: float = 0.0, maximum: float = 1.0) -> float:
    try:
        value = float(values[key])
    except ValueError as exc:
        raise ConfError(f"{key} must be numeric") from exc
    if not minimum <= value <= maximum:
        raise ConfError(f"{key} must be between {minimum} and {maximum}")
    return value


def _positive_float(values: dict[str, str], key: str) -> float:
    try:
        value = float(values[key])
    except ValueError as exc:
        raise ConfError(f"{key} must be numeric") from exc
    if value <= 0:
        raise ConfError(f"{key} must be > 0")
    return value


def _nonnegative_float(values: dict[str, str], key: str) -> float:
    try:
        value = float(values[key])
    except ValueError as exc:
        raise ConfError(f"{key} must be numeric") from exc
    if value < 0:
        raise ConfError(f"{key} must be >= 0")
    return value


def _path(value: str, base: Path) -> str:
    if not value:
        return ""
    path = Path(value).expanduser()
    return str((path if path.is_absolute() else base / path).resolve())


def project_from_conf(path: str | Path) -> tuple[dict[str, Any], str]:
    values = read_conf(path)
    base = Path(values["_CONFIG_PATH"]).parent
    for key, accepted in {
        "LIBRARY_LAYOUT": {"single_end", "se", "paired_end", "pe"},
        "MAPPING_POLICY": {"unique_primary"},
        "END_SOFT_CLIP_POLICY": {"exclude_and_report"},
        "PAIRING_MODE": {"auto", "none", "required"},
        "INCOMPLETE_PAIR_ACTION": {"error", "unpaired"},
        "AMBIGUOUS_GENE_POLICY": {"exclude_statistics"},
        "PAS_MASK_RESCUE_TIER": {"core", "core_plus_rescue"},
        "PAS_DISCOVERY_THRESHOLD_OPERATOR": {"greater_than"},
        "APA_B_ENDPOINT_SOURCE": {"auto", "exact_ends", "bam"},
        "FASTQ_SCREEN_MISSING_ACTION": {"warn", "error"},
    }.items():
        if values[key].lower() not in accepted:
            raise ConfError(f"Unsupported {key}: {values[key]}")
    if _int(values, "PAS_DISCOVERY_WINDOW_NT") != 30 or _int(values, "PAS_DISCOVERY_ROUNDS") != 2:
        raise ConfError("Faithful Mcell2019 mode requires 30-nt windows and exactly two discovery rounds")
    if _bool(values, "UMI_PRESENT"):
        raise ConfError("UMI_PRESENT must be false")
    if _bool(values, "RUN_APA_B") and not _bool(values, "APA_B_PILOT_ACCEPTED"):
        raise ConfError("RUN_APA_B=true requires explicit APA_B_PILOT_ACCEPTED=true")
    if _bool(values, "RUN_APA_B") and not values["APA_B_COMMAND_TEMPLATE"].strip():
        raise ConfError("RUN_APA_B=true requires APA_B_COMMAND_TEMPLATE")
    if _bool(values, "RUN_APA_B") and not values["APA_B_VALIDATION_MANIFEST"].strip():
        raise ConfError("RUN_APA_B=true requires APA_B_VALIDATION_MANIFEST")
    if _bool(values, "RUN_APA_B") and not values["APA_B_INSTALLATION_MANIFEST"].strip():
        raise ConfError("RUN_APA_B=true requires APA_B_INSTALLATION_MANIFEST")
    apa_b_threads = _int(values, "APA_B_THREADS")
    for key in ("APA_B_ENDPOINT_PARALLEL_JOBS", "APA_B_CLUSTER_PARALLEL_JOBS", "APA_B_DEEPIP_THREADS"):
        value = _int(values, key)
        if value < 1 or value > apa_b_threads:
            raise ConfError(f"{key} must be between 1 and APA_B_THREADS ({apa_b_threads})")
    if not (_bool(values, "ENRICHMENT_ORA") or _bool(values, "ENRICHMENT_GSEA")):
        raise ConfError("At least one of ENRICHMENT_ORA or ENRICHMENT_GSEA must be true")
    if _bool(values, "RUN_RSEQC") and not any(_bool(values, key) for key in (
        "RSEQC_INFER_EXPERIMENT", "RSEQC_READ_DISTRIBUTION", "RSEQC_GENE_BODY_COVERAGE",
    )):
        raise ConfError("RUN_RSEQC=true requires at least one enabled RSeQC analysis")
    if _int(values, "RSEQC_MIN_TRANSCRIPT_LENGTH") < 100:
        raise ConfError("RSEQC_MIN_TRANSCRIPT_LENGTH must be at least 100")
    minimum_geneset = _int(values, "ENRICHMENT_MIN_GENESET_SIZE")
    maximum_geneset = _int(values, "ENRICHMENT_MAX_GENESET_SIZE")
    if minimum_geneset > maximum_geneset:
        raise ConfError("ENRICHMENT_MIN_GENESET_SIZE must not exceed ENRICHMENT_MAX_GENESET_SIZE")
    if (_bool(values, "RUN_TRACKS") and not _bool(values, "GENERATE_BIGWIGS")
            and not (_bool(values, "RETAIN_BEDGRAPH") or _bool(values, "KEEP_BEDGRAPHS"))):
        raise ConfError("Track generation requires BigWig and/or retained bedGraph output")
    ucsc_url = values["UCSC_BIGDATA_URL_PREFIX"].rstrip("/")
    if any(character in ucsc_url for character in "[]()"):
        raise ConfError("UCSC_BIGDATA_URL_PREFIX must be a plain URL, not Markdown link syntax")
    if ucsc_url and not re.match(r"^https?://[^\s]+$", ucsc_url):
        raise ConfError("UCSC_BIGDATA_URL_PREFIX must be empty or an http(s) URL")
    if not re.fullmatch(r"-?[0-9]+(?:\.[0-9]+)?:-?[0-9]+(?:\.[0-9]+)?", values["UCSC_VIEW_LIMITS"]):
        raise ConfError("UCSC_VIEW_LIMITS must have numeric min:max syntax")
    condition_order = [item.strip() for item in values["CONDITION_ORDER"].split(",") if item.strip()]
    if len(condition_order) != len(set(condition_order)):
        raise ConfError("CONDITION_ORDER must not contain duplicate conditions")
    family_keys = ("GENERATE_ALL_READ_TRACKS", "GENERATE_EXACT_END_TRACKS", "GENERATE_FILTERED_END_TRACKS",
                   "GENERATE_REJECTED_END_TRACKS", "GENERATE_ACTIVE_PAS_TRACKS")
    normalization_keys = ("GENERATE_RAW_TRACKS", "GENERATE_CPM_TRACKS", "GENERATE_DESEQ2_FINAL_TRACKS",
                          "GENERATE_DESEQ2_ROBUST_CPM_FINAL_TRACKS")
    if _bool(values, "RUN_TRACKS") and not any(_bool(values, key) for key in family_keys):
        raise ConfError("RUN_TRACKS=true requires at least one enabled track family")
    if _bool(values, "RUN_TRACKS") and not any(_bool(values, key) for key in normalization_keys):
        raise ConfError("RUN_TRACKS=true requires at least one enabled track normalization")
    nonfinal_family = any(_bool(values, key) for key in
                          ("GENERATE_ALL_READ_TRACKS", "GENERATE_EXACT_END_TRACKS", "GENERATE_REJECTED_END_TRACKS"))
    base_norm = _bool(values, "GENERATE_RAW_TRACKS") or _bool(values, "GENERATE_CPM_TRACKS")
    final_family = _bool(values, "GENERATE_FILTERED_END_TRACKS") or _bool(values, "GENERATE_ACTIVE_PAS_TRACKS")
    final_norm = (_bool(values, "GENERATE_DESEQ2_FINAL_TRACKS") or
                  _bool(values, "GENERATE_DESEQ2_ROBUST_CPM_FINAL_TRACKS"))
    if _bool(values, "RUN_TRACKS") and nonfinal_family and not base_norm:
        raise ConfError("C0/C1/C2R track families require raw and/or CPM normalization")
    if _bool(values, "RUN_TRACKS") and final_norm and not final_family:
        raise ConfError("DESeq2 final-track normalizations require filtered-end and/or active-PAS tracks")

    references = {
        "GRCh38": {
            "species": "human", "assembly": "GRCh38", "release": "GENCODE_v42",
            "star_index": _path(values.get("HG38_STAR_INDEX", ""), base),
            "fasta": _path(values.get("HG38_FASTA", ""), base),
            "gtf": _path(values.get("HG38_GTF", ""), base),
            "chrom_sizes": _path(values.get("HG38_CHROM_SIZES", ""), base),
            "pas_atlas": _path(values.get("HG38_PAS_ATLAS", ""), base),
            "rseqc_bed": _path(values.get("HG38_RSEQC_BED", ""), base),
        },
        "GRCm39": {
            "species": "mouse", "assembly": "GRCm39", "release": "GENCODE_vM31",
            "star_index": _path(values.get("MM39_STAR_INDEX", ""), base),
            "fasta": _path(values.get("MM39_FASTA", ""), base),
            "gtf": _path(values.get("MM39_GTF", ""), base),
            "chrom_sizes": _path(values.get("MM39_CHROM_SIZES", ""), base),
            "pas_atlas": _path(values.get("MM39_PAS_ATLAS", ""), base),
            "rseqc_bed": _path(values.get("MM39_RSEQC_BED", ""), base),
        },
    }
    project: dict[str, Any] = {
        "_config_path": values["_CONFIG_PATH"],
        "_samplesheet_path": _path(values["SAMPLESHEET"], base),
        "_config_format": "conf",
        "project_id": values["PROJECT_ID"],
        "output_dir": _path(values["OUTPUT_DIR"], base),
        "tmp_dir": _path(values["TMP_DIR"], base),
        "condition_order": condition_order,
        "design": values["DESIGN"],
        "protocol": {
            "profile": values["LIBRARY_PROTOCOL"].lower(), "has_umi": False,
            "library_layout": "PE" if values["LIBRARY_LAYOUT"].lower() in {"paired_end", "pe"} else "SE",
            "end_defining_mate": "R1",
            "retain_duplicate_flagged_reads": True, "mapping_policy": values["MAPPING_POLICY"],
            "end_soft_clip_policy": values["END_SOFT_CLIP_POLICY"],
            "orientation_min_fraction": _float(values, "ORIENTATION_MIN_FRACTION"),
        },
        "preprocessing": {
            "bbduk_reference": (
                values["BBDUK_REFERENCE"] if values["BBDUK_REFERENCE"] == "adapters,polyA_T"
                else _path(values["BBDUK_REFERENCE"], base)
            ),
            "trim_quality": _int(values, "TRIM_QUALITY", 0),
            "minimum_length": _int(values, "MINIMUM_READ_LENGTH"),
            "pe_r2_trim_5p": _int(values, "PE_R2_TRIM_5P", 0),
            "fastq_screen": {
                "enabled": _bool(values, "RUN_FASTQ_SCREEN"),
                "config": _path(values["FASTQ_SCREEN_CONFIG"], base),
                "missing_action": values["FASTQ_SCREEN_MISSING_ACTION"].lower(),
                "subset": _int(values, "FASTQ_SCREEN_SUBSET"),
            },
        },
        "modules": {
            "gene_expression": _bool(values, "RUN_GENE_EXPRESSION"),
            "apa_a": _bool(values, "RUN_APA_A_MCELL2019"),
            "apa_b": _bool(values, "RUN_APA_B"), "tracks": _bool(values, "RUN_TRACKS"),
            "rseqc": _bool(values, "RUN_RSEQC"),
            "dge_enrichment": _bool(values, "RUN_DGE_ENRICHMENT"),
            "apa_enrichment": _bool(values, "RUN_APA_ENRICHMENT"),
        },
        "references": references,
        "statistics": {
            "min_replicates": _int(values, "MIN_REPLICATES_PER_CONDITION", 2),
            "ambiguous_gene_policy": values["AMBIGUOUS_GENE_POLICY"],
            "pairing": {
                "mode": values["PAIRING_MODE"], "subject_column": values["PAIRING_COLUMN"],
                "paired_design": f"~ {values['PAIRING_COLUMN']} + condition",
                "incomplete_pair_action": values["INCOMPLETE_PAIR_ACTION"],
            },
        },
        "apa_a": {
            "method": "mcell2019", "mapping_policy": values["MAPPING_POLICY"],
            "soft_clip_policy": values["END_SOFT_CLIP_POLICY"],
            "internal_priming_consecutive_bases": _int(values, "INTERNAL_PRIMING_CONSECUTIVE_BASES"),
            "internal_priming_window_nt": _int(values, "INTERNAL_PRIMING_WINDOW_NT"),
            "internal_priming_min_bases_in_window": _int(values, "INTERNAL_PRIMING_MIN_BASES_IN_WINDOW"),
            "mask_rescue_tier": values["PAS_MASK_RESCUE_TIER"],
            "discovery_window_nt": _int(values, "PAS_DISCOVERY_WINDOW_NT"),
            "discovery_threshold": _positive_float(values, "PAS_DISCOVERY_THRESHOLD"),
            "discovery_threshold_operator": values["PAS_DISCOVERY_THRESHOLD_OPERATOR"],
            "discovery_rounds": _int(values, "PAS_DISCOVERY_ROUNDS"),
            "gene_downstream_extension_nt": _int(values, "GENE_DOWNSTREAM_EXTENSION_NT", 0),
        },
        "apa_b": {
            "enabled": _bool(values, "RUN_APA_B"), "pilot_accepted": _bool(values, "APA_B_PILOT_ACCEPTED"),
            "command_template": values["APA_B_COMMAND_TEMPLATE"],
            "validation_manifest": _path(values["APA_B_VALIDATION_MANIFEST"], base),
            "installation_manifest": _path(values["APA_B_INSTALLATION_MANIFEST"], base),
            "endpoint_source": values["APA_B_ENDPOINT_SOURCE"].lower(),
        },
        "reporting": {"fdr": _float(values, "FDR"), "min_abs_delta_pau": _float(values, "MIN_ABS_DELTA_PAU")},
        "rseqc": {
            "enabled": _bool(values, "RUN_RSEQC"),
            "infer_experiment": _bool(values, "RSEQC_INFER_EXPERIMENT"),
            "read_distribution": _bool(values, "RSEQC_READ_DISTRIBUTION"),
            "gene_body_coverage": _bool(values, "RSEQC_GENE_BODY_COVERAGE"),
            "multiqc": _bool(values, "RSEQC_MULTIQC"),
            "sample_reads": _int(values, "RSEQC_SAMPLE_READS"),
            "minimum_transcript_length": _int(values, "RSEQC_MIN_TRANSCRIPT_LENGTH"),
        },
        "enrichment": {
            "ora": _bool(values, "ENRICHMENT_ORA"), "gsea": _bool(values, "ENRICHMENT_GSEA"),
            "go": _bool(values, "ENRICHMENT_GO"), "reactome": _bool(values, "ENRICHMENT_REACTOME"),
            "hallmarks": _bool(values, "ENRICHMENT_HALLMARKS"),
            "kegg": _bool(values, "ENRICHMENT_KEGG"),
            "rich_plots": _bool(values, "ENRICHMENT_RICH_PLOTS"),
            "network_max_terms": _int(values, "ENRICHMENT_NETWORK_MAX_TERMS"),
            "network_max_genes": _int(values, "ENRICHMENT_NETWORK_MAX_GENES"),
            "padj": _float(values, "ENRICHMENT_PADJ"),
            "dge_min_abs_lfc": _nonnegative_float(values, "ENRICHMENT_DGE_MIN_ABS_LFC"),
            "apa_min_abs_delta_pau": _nonnegative_float(values, "ENRICHMENT_APA_MIN_ABS_DELTA_PAU"),
            "min_geneset_size": minimum_geneset, "max_geneset_size": maximum_geneset,
        },
        "tracks": {
            "early_c0": _bool(values, "GENERATE_EARLY_C0_TRACKS"),
            "families": {
                "all_reads": _bool(values, "GENERATE_ALL_READ_TRACKS"),
                "exact_ends": _bool(values, "GENERATE_EXACT_END_TRACKS"),
                "filtered_ends": _bool(values, "GENERATE_FILTERED_END_TRACKS"),
                "rejected_ends": _bool(values, "GENERATE_REJECTED_END_TRACKS"),
                "active_pas": _bool(values, "GENERATE_ACTIVE_PAS_TRACKS"),
            },
            "normalizations": {
                "raw": _bool(values, "GENERATE_RAW_TRACKS"), "cpm": _bool(values, "GENERATE_CPM_TRACKS"),
                "deseq2": _bool(values, "GENERATE_DESEQ2_FINAL_TRACKS"),
                "robust_cpm": _bool(values, "GENERATE_DESEQ2_ROBUST_CPM_FINAL_TRACKS"),
            },
            "generate_bigwigs": _bool(values, "GENERATE_BIGWIGS"),
            "retain_bedgraph": _bool(values, "RETAIN_BEDGRAPH") or _bool(values, "KEEP_BEDGRAPHS"),
            "ucsc_bigdata_url_prefix": ucsc_url,
            "ucsc_negate_minus_tracks": _bool(values, "UCSC_NEGATE_MINUS_TRACKS"),
            "ucsc_view_limits": values["UCSC_VIEW_LIMITS"],
        },
        "cleanup": {
            "enabled": _bool(values, "CLEANUP_INTERMEDIATES"),
            "keep_trimmed_fastq": _bool(values, "KEEP_TRIMMED_FASTQ"),
            "keep_lane_bams": _bool(values, "KEEP_LANE_BAMS"),
            "keep_apa_sample_extraction": _bool(values, "KEEP_APA_SAMPLE_EXTRACTION"),
            "keep_track_strand_bams": _bool(values, "KEEP_TRACK_STRAND_BAMS"),
            "keep_track_bedgraphs": _bool(values, "KEEP_BEDGRAPHS") or _bool(values, "RETAIN_BEDGRAPH"),
        },
        "resources": {
            "total_threads": _int(values, "MAX_TOTAL_THREADS"), "total_memory_gb": _int(values, "MAX_TOTAL_MEMORY_GB"),
            "temporary_directory": _path(values["TMP_DIR"], base),
            "downstream": {
                "parallel_modules": (
                    _int(values, "DOWNSTREAM_MODULE_PARALLEL_JOBS")
                    if _bool(values, "PARALLEL_DOWNSTREAM_MODULES") else 1
                ),
            },
            "preprocess": {
                "trim_parallel_jobs": _int(values, "PREPROCESS_PARALLEL_JOBS"),
                "star_parallel_jobs": _int(values, "STAR_PARALLEL_JOBS"),
                "fastqc_threads": _int(values, "FASTQC_THREADS"),
                "bbduk_threads": _int(values, "BBDUK_THREADS"), "bbduk_memory_gb": _int(values, "BBDUK_MEMORY_GB"),
                "star_threads": _int(values, "STAR_THREADS"), "star_memory_gb": _int(values, "STAR_MEMORY_GB"),
                "samtools_threads": _int(values, "SAMTOOLS_THREADS"),
                "samtools_sort_memory_per_thread_gb": _int(values, "SAMTOOLS_SORT_MEMORY_PER_THREAD_GB"),
                "merge_parallel_jobs": _int(values, "SAMPLE_MERGE_PARALLEL_JOBS"), "merge_memory_gb": 16,
                "fastq_screen_parallel_jobs": _int(values, "FASTQ_SCREEN_PARALLEL_JOBS"),
                "fastq_screen_threads": _int(values, "FASTQ_SCREEN_THREADS"),
                "fastq_screen_memory_gb": _int(values, "FASTQ_SCREEN_MEMORY_GB"),
            },
            "rseqc": {
                "parallel_jobs": _int(values, "RSEQC_PARALLEL_JOBS"),
                "threads": 1,
                "memory_gb": _int(values, "RSEQC_MEMORY_GB"),
            },
            "dge": {"featurecounts_threads": _int(values, "SAMTOOLS_THREADS"), "featurecounts_memory_gb": 16,
                    "contrast_parallel_jobs": _int(values, "DGE_CONTRAST_PARALLEL_JOBS"), "contrast_threads": 1, "contrast_memory_gb": 16},
            "apa_a": {"extraction_parallel_jobs": _int(values, "END_EXTRACTION_PARALLEL_JOBS"), "extraction_threads": 1,
                      "extraction_memory_gb": 4,
                      "contrast_parallel_jobs": _int_or(
                          values, "APA_A_CONTRAST_PARALLEL_JOBS", "APA_CONTRAST_PARALLEL_JOBS"),
                      "contrast_threads": 1, "contrast_memory_gb": 16},
            "apa_b": {"engine_threads": _int(values, "APA_B_THREADS"), "engine_memory_gb": 24,
                      "endpoint_parallel_jobs": _int(values, "APA_B_ENDPOINT_PARALLEL_JOBS"),
                      "cluster_parallel_jobs": _int(values, "APA_B_CLUSTER_PARALLEL_JOBS"),
                      "deepip_threads": _int(values, "APA_B_DEEPIP_THREADS"), "sample_memory_gb": 4,
                      "contrast_parallel_jobs": _int_or(
                          values, "APA_B_CONTRAST_PARALLEL_JOBS", "APA_CONTRAST_PARALLEL_JOBS"),
                      "contrast_threads": 1,
                      "contrast_memory_gb": 16},
            "enrichment": {"parallel_jobs": _int(values, "ENRICHMENT_PARALLEL_JOBS"),
                           "threads": 1, "memory_gb": 16},
            "tracks": {"parallel_jobs": _int(values, "TRACK_PARALLEL_JOBS"),
                       "samtools_threads": _int(values, "TRACK_THREADS"), "memory_gb": 8},
        },
    }
    fixed_method = {
        "INTERNAL_PRIMING_CONSECUTIVE_BASES": "6", "INTERNAL_PRIMING_WINDOW_NT": "10",
        "INTERNAL_PRIMING_MIN_BASES_IN_WINDOW": "7", "PAS_DISCOVERY_WINDOW_NT": "30",
        "PAS_DISCOVERY_THRESHOLD": "30", "PAS_DISCOVERY_THRESHOLD_OPERATOR": "greater_than",
        "PAS_DISCOVERY_ROUNDS": "2", "GENE_DOWNSTREAM_EXTENSION_NT": "6000",
    }
    project["apa_a"]["method_status"] = (
        "Mcell2019_faithful" if all(values[key] == expected for key, expected in fixed_method.items())
        else "Mcell2019_modified"
    )
    return project, _path(values["SAMPLESHEET"], base)
