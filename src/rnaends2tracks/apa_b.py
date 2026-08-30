from __future__ import annotations

import csv
import json
import re
import shlex
from pathlib import Path

from .config import RunPlan, signature_for
from .external import event, require_tools, run
from .mcell2019 import load_gene_models
from .receipts import receipt_valid, write_receipt
from .statistics import run_r_contrasts

REQUIRED_CATALOG = {"pas_id", "gene_id", "chrom", "start", "end", "strand", "feature_class"}
REQUIRED_NA_AUDIT = {
    "status", "screening_tests", "screening_na", "screening_na_fraction",
    "confirmation_tests", "confirmation_na", "confirmation_na_fraction",
    "excluded_genes_fewer_than_two_testable_sites", "stageR_input_genes",
    "stageR_input_sites", "stageR_adjusted_na", "stageR_adjusted_na_fraction",
    "na_policy",
}
REQUIRED_FIT_AUDIT = {
    "contrast_id", "status", "fit_policy", "multifactor", "one_way",
    "random_seed", "add_uniform_used", "primary_error",
}


def _validation_manifest(
    path: Path, genomes: list[str], library_protocol: str = "quantseq_rev_v2_se",
) -> dict[str, object]:
    if not path.is_file():
        raise RuntimeError(f"APA-B validation manifest is unavailable: {path}")
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError) as exc:
        raise RuntimeError(f"APA-B validation manifest is invalid JSON: {path}") from exc
    if payload.get("schema_version") != 1 or payload.get("status") != "accepted":
        raise RuntimeError("APA-B validation manifest must have schema_version=1 and status=accepted")
    engine = payload.get("engine", {})
    workflow_adapter = payload.get("workflow_adapter", {})
    model = payload.get("model", {})
    models = payload.get("models", {})
    environment = payload.get("environment", {})
    if not isinstance(engine, dict) or not re.fullmatch(r"[0-9a-fA-F]{7,40}", str(engine.get("source_commit", ""))):
        raise RuntimeError("APA-B validation manifest requires a pinned engine.source_commit")
    if (not isinstance(workflow_adapter, dict)
            or not re.fullmatch(r"[0-9a-fA-F]{7,40}", str(workflow_adapter.get("source_commit", "")))):
        raise RuntimeError("APA-B validation manifest requires a pinned workflow_adapter.source_commit")
    if not model and isinstance(models, dict) and models:
        model = next(iter(models.values()))
    for label, record in (("model/models", model), ("environment", environment)):
        if not isinstance(record, dict) or not re.fullmatch(r"[0-9a-fA-F]{64}", str(record.get("sha256", ""))):
            raise RuntimeError(f"APA-B validation manifest requires {label}.sha256")
    if isinstance(models, dict) and models:
        assembly_species = {"GRCh38": "human", "GRCm39": "mouse"}
        for genome in genomes:
            species = assembly_species.get(genome)
            record = models.get(species, {})
            if not isinstance(record, dict) or not re.fullmatch(r"[0-9a-fA-F]{64}", str(record.get("sha256", ""))):
                raise RuntimeError(f"APA-B validation manifest lacks a pinned {species} model for {genome}")
    if payload.get("umi_present") is not False or payload.get("coordinate_deduplication") is not False:
        raise RuntimeError("APA-B validation must explicitly record no UMI and no coordinate deduplication")
    if payload.get("quantseq_rev_adaptation") != "genomewide_no_tail_weighted_PAC":
        raise RuntimeError("APA-B validation manifest covers a different QuantSeq REV adaptation")
    protocols = payload.get("library_protocols")
    if protocols is None:
        # Accepted manifests created before paired-end support represented the
        # validated QuantSeq REV V2 SE contract implicitly. Keep that narrow
        # compatibility path; PE always requires a newly reviewed PE canary.
        protocols = ["quantseq_rev_v2_se"]
        payload["library_protocols"] = protocols
    if not isinstance(protocols, list) or library_protocol not in protocols:
        raise RuntimeError(f"APA-B validation does not cover {library_protocol}")
    missing = sorted(set(genomes).difference(map(str, payload.get("assemblies", []))))
    if missing:
        raise RuntimeError("APA-B validation does not cover assemblies: " + ", ".join(missing))
    pilot = payload.get("pilot", {})
    if (not isinstance(pilot, dict) or pilot.get("synthetic_pass") is not True
            or pilot.get("c1_c1s_reuse_equivalent") is not True):
        raise RuntimeError("APA-B validation requires passing synthetic and C1+C1S equivalence pilots")
    real_canaries = pilot.get("real_quantseq_rev_canaries", {})
    if not isinstance(real_canaries, dict) or any(real_canaries.get(genome) != "PASS" for genome in genomes):
        raise RuntimeError("APA-B validation requires a passing real QuantSeq REV canary for every assembly")
    for field in ("reviewed_by", "accepted_at"):
        if not str(payload.get(field, "")).strip():
            raise RuntimeError(f"APA-B validation manifest requires {field}")
    return payload


def _validate_engine_provenance(path: Path, accepted: dict[str, object], assembly: str) -> None:
    if not path.is_file():
        raise RuntimeError(f"APA-B adapter did not create engine provenance: {path}")
    observed = json.loads(path.read_text(encoding="utf-8"))
    if observed.get("assembly") != assembly:
        raise RuntimeError(f"APA-B engine provenance assembly mismatch: {observed.get('assembly')} != {assembly}")
    for section, field in (("workflow_adapter", "source_commit"), ("engine", "source_commit"),
                           ("model", "sha256"), ("environment", "sha256")):
        if section == "model" and isinstance(accepted.get("models"), dict):
            expected_section = accepted["models"].get(str(observed.get("species", "")), accepted.get("model", {}))
        else:
            expected_section = accepted.get(section, {})
        observed_section = observed.get(section, {})
        if not isinstance(expected_section, dict) or not isinstance(observed_section, dict):
            raise RuntimeError(f"APA-B engine provenance lacks {section}")
        if str(observed_section.get(field, "")).lower() != str(expected_section.get(field, "")).lower():
            raise RuntimeError(f"APA-B engine provenance does not match accepted {section}.{field}")
    if observed.get("umi_present") is not False or observed.get("coordinate_deduplication") is not False:
        raise RuntimeError("APA-B engine provenance violates the accepted no-UMI/no-dedup contract")
    if observed.get("library_protocol") not in accepted.get("library_protocols", []):
        raise RuntimeError("APA-B engine provenance reports a library protocol outside the accepted pilot")


def _header(path: Path) -> list[str]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return next(csv.reader(handle, delimiter="\t"), [])


def _validate_polyaseqtrap_outputs(catalog: Path, counts: Path, deepip: Path, sample_ids: list[str]) -> None:
    for path in (catalog, counts, deepip):
        if not path.is_file() or path.stat().st_size == 0:
            raise RuntimeError(f"PolyAseqTrap adapter did not create required output: {path}")
    catalog_header = set(_header(catalog))
    if not REQUIRED_CATALOG.issubset(catalog_header):
        raise RuntimeError("APA-B catalog requires columns: " + ", ".join(sorted(REQUIRED_CATALOG)))
    counts_header = _header(counts)
    expected = ["pas_id", *sample_ids]
    if counts_header != expected:
        raise RuntimeError(f"APA-B count columns must exactly equal {expected}; observed {counts_header}")
    for path in (catalog, counts):
        with path.open(encoding="utf-8", errors="ignore") as handle:
            if "APA_A_" in handle.read(10000):
                raise RuntimeError("APA-B output contains APA-A identifiers; independent catalogs are mandatory")
    with counts.open(encoding="utf-8", newline="") as handle:
        for line_number, row in enumerate(csv.DictReader(handle, delimiter="\t"), start=2):
            for sample_id in sample_ids:
                try:
                    value = int(row[sample_id])
                except ValueError as exc:
                    raise RuntimeError(f"Non-integer APA-B count at line {line_number}, sample {sample_id}") from exc
                if value < 0:
                    raise RuntimeError(f"Negative APA-B count at line {line_number}, sample {sample_id}")


def _validate_na_audit(path: Path) -> None:
    if not path.is_file() or path.stat().st_size == 0:
        raise RuntimeError(f"APA-B contrast NA audit is unavailable: {path}")
    with path.open(encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if len(rows) != 1 or not REQUIRED_NA_AUDIT.issubset(rows[0]):
        raise RuntimeError(f"APA-B contrast NA audit has an invalid schema: {path}")
    row = rows[0]
    if row["status"] not in {"PASS", "WARN_UNTESTABLE_PVALUES"}:
        raise RuntimeError(f"APA-B contrast NA audit has an invalid status: {path}")
    integer_fields = (
        "screening_tests", "screening_na", "confirmation_tests", "confirmation_na",
        "excluded_genes_fewer_than_two_testable_sites", "stageR_input_genes",
        "stageR_input_sites", "stageR_adjusted_na",
    )
    fraction_fields = (
        "screening_na_fraction", "confirmation_na_fraction", "stageR_adjusted_na_fraction",
    )
    try:
        integers = {field: int(row[field]) for field in integer_fields}
        fractions = {field: float(row[field]) for field in fraction_fields}
    except ValueError as exc:
        raise RuntimeError(f"APA-B contrast NA audit contains a non-numeric value: {path}") from exc
    if any(value < 0 for value in integers.values()):
        raise RuntimeError(f"APA-B contrast NA audit contains a negative count: {path}")
    if any(not 0 <= value <= 1 for value in fractions.values()):
        raise RuntimeError(f"APA-B contrast NA audit contains an invalid fraction: {path}")
    if (integers["screening_na"] > integers["screening_tests"]
            or integers["confirmation_na"] > integers["confirmation_tests"]
            or integers["stageR_input_genes"] > integers["screening_tests"]
            or integers["stageR_input_sites"] > integers["confirmation_tests"]
            or integers["stageR_adjusted_na"] > integers["confirmation_tests"]):
        raise RuntimeError(f"APA-B contrast NA audit contains inconsistent counts: {path}")
    if integers["screening_tests"] == 0 or integers["confirmation_tests"] == 0:
        raise RuntimeError(f"APA-B contrast NA audit reports no DRIMSeq tests: {path}")
    expected_fractions = {
        "screening_na_fraction": integers["screening_na"] / integers["screening_tests"],
        "confirmation_na_fraction": integers["confirmation_na"] / integers["confirmation_tests"],
        "stageR_adjusted_na_fraction": integers["stageR_adjusted_na"] / integers["confirmation_tests"],
    }
    if any(abs(fractions[field] - expected) > 1e-9 for field, expected in expected_fractions.items()):
        raise RuntimeError(f"APA-B contrast NA audit contains inconsistent fractions: {path}")
    if integers["stageR_input_genes"] == 0 or integers["stageR_input_sites"] == 0:
        raise RuntimeError(f"APA-B contrast NA audit reports no testable stageR hypotheses: {path}")
    if row["na_policy"] != "untestable hypotheses remain NA and cannot be significant":
        raise RuntimeError(f"APA-B contrast NA audit reports an unsupported NA policy: {path}")


def _validate_fit_audit(path: Path) -> None:
    if not path.is_file() or path.stat().st_size == 0:
        raise RuntimeError(f"APA-B contrast fit audit is unavailable: {path}")
    with path.open(encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if len(rows) != 1 or not REQUIRED_FIT_AUDIT.issubset(rows[0]):
        raise RuntimeError(f"APA-B contrast fit audit has an invalid schema: {path}")
    row = rows[0]
    policy = row["fit_policy"]
    expected_status = {
        "standard": "PASS",
        "deterministic_add_uniform_retry": "WARN_NUMERIC_RETRY",
    }
    if policy not in expected_status or row["status"] != expected_status[policy]:
        raise RuntimeError(f"APA-B contrast fit audit has an invalid policy/status: {path}")
    booleans = {"TRUE": True, "FALSE": False}
    try:
        multifactor = booleans[row["multifactor"]]
        one_way = booleans[row["one_way"]]
        add_uniform = booleans[row["add_uniform_used"]]
        seed = int(row["random_seed"])
    except (KeyError, ValueError) as exc:
        raise RuntimeError(f"APA-B contrast fit audit has invalid typed values: {path}") from exc
    if seed <= 0 or one_way == multifactor:
        raise RuntimeError(f"APA-B contrast fit audit has inconsistent design metadata: {path}")
    if policy == "standard" and (add_uniform or row["primary_error"]):
        raise RuntimeError(f"APA-B standard fit audit unexpectedly records a retry: {path}")
    if policy == "deterministic_add_uniform_retry" and (
        not multifactor or not add_uniform or not row["primary_error"]
    ):
        raise RuntimeError(f"APA-B numerical-retry audit is incomplete: {path}")


def apa_b(plan: RunPlan, results: Path, script_root: Path, dry_run: bool = False, force: bool = False) -> None:
    settings = plan.project.get("apa_b", {})
    module_dir = results / "07_apa_b"
    log_dir = results / "logs"
    module_dir.mkdir(parents=True, exist_ok=True)
    if settings.get("enabled", False) is False:
        event(log_dir, "apa_b", "disabled", "APA-B disabled in project configuration")
        return
    if settings.get("pilot_accepted") is not True:
        raise RuntimeError(
            "APA-B is pilot-gated. Set apa_b.pilot_accepted: true only after exact-end, no-dedup, "
            "DeepIP, and count-conservation acceptance criteria pass."
        )
    template = str(settings.get("command_template", "")).strip()
    if not template:
        raise RuntimeError("apa_b.command_template is required for the pinned local PolyAseqTrap installation")
    validation_path = Path(str(settings.get("validation_manifest", "")))
    installation_path = Path(str(settings.get("installation_manifest", "")))
    library_protocol = str(plan.project.get("protocol", {}).get("profile", "quantseq_rev_v2_se"))
    accepted_validation = _validation_manifest(validation_path, list(plan.references), library_protocol)
    bams = [results / "02_alignment" / sample["sample_id"] / f"{sample['sample_id']}.bam" for sample in plan.samples]
    endpoint_inputs = [
        results / "03_exact_ends" / sample["genome"] / sample["sample_id"] / filename
        for sample in plan.samples
        for filename in ("C1_exact_ends.tsv.gz", "C1S_uncertain_ends.tsv.gz", "end_audit.json")
    ]
    endpoint_receipts = [
        results / "03_exact_ends" / sample["genome"] / sample["sample_id"] / ".receipt" / "run_receipt.json"
        for sample in plan.samples
    ]
    ref_inputs = [Path(plan.reference_for(genome)[key]) for genome in plan.references for key in ("fasta", "gtf")]
    reusable_inputs = [*endpoint_inputs, *endpoint_receipts]
    signature_inputs = [*ref_inputs, validation_path, installation_path]
    if settings.get("endpoint_source", "auto") != "bam" and all(path.is_file() for path in reusable_inputs):
        signature_inputs.extend(reusable_inputs)
    else:
        signature_inputs.extend(bams)
    signature = signature_for(signature_inputs, {
        "module": "apa_b", "settings": settings, "samples": plan.samples,
        "contrasts": plan.contrasts, "reporting": plan.project.get("reporting", {}),
    }) if not dry_run else "dry-run"
    if not force and not dry_run and receipt_valid(module_dir, signature):
        event(log_dir, "apa_b", "skipped", "Valid matching receipt")
        return
    outputs: list[Path] = []
    if not dry_run:
        require_tools(["Rscript"])
    for genome in plan.references:
        samples = [sample for sample in plan.samples if sample["genome"] == genome]
        contrasts = [contrast for contrast in plan.contrasts if contrast["genome"] == genome]
        if not samples or not contrasts:
            continue
        reference = plan.reference_for(genome)
        genome_dir = module_dir / genome
        genome_dir.mkdir(parents=True, exist_ok=True)
        manifest = genome_dir / "bam_manifest.tsv"
        genome_bams = [results / "02_alignment" / sample["sample_id"] / f"{sample['sample_id']}.bam" for sample in samples]
        with manifest.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["sample_id", "bam", "library_layout", "end_defining_mate",
                             "c1", "c1s", "exact_end_audit", "exact_end_receipt"])
            for sample, bam in zip(samples, genome_bams):
                exact_root = results / "03_exact_ends" / genome / sample["sample_id"]
                writer.writerow([
                    sample["sample_id"], bam, sample.get("library_layout", "SE"), "R1",
                    exact_root / "C1_exact_ends.tsv.gz",
                    exact_root / "C1S_uncertain_ends.tsv.gz",
                    exact_root / "end_audit.json",
                    exact_root / ".receipt" / "run_receipt.json",
                ])
        replacements = {"bam_manifest": str(manifest), "fasta": reference["fasta"],
                        "gtf": reference["gtf"], "outdir": str(genome_dir),
                        "species": reference["species"], "assembly": reference["assembly"],
                        "threads": str(plan.project["resources"]["apa_b"]["engine_threads"]),
                        "endpoint_source": str(settings.get("endpoint_source", "auto")),
                        "endpoint_workers": str(plan.project["resources"]["apa_b"]["endpoint_parallel_jobs"]),
                        "cluster_workers": str(plan.project["resources"]["apa_b"]["cluster_parallel_jobs"]),
                        "deepip_threads": str(plan.project["resources"]["apa_b"]["deepip_threads"]),
                        "library_protocol": library_protocol,
                        "validation_manifest": str(validation_path),
                        "installation_manifest": str(settings.get("installation_manifest", "")),
                        "apa_b_executable": str(
                            Path(str(settings.get("installation_manifest", ""))).parent
                            / "bin" / "rna-ends2tracks-apa-b"
                        )}
        if template.lower() == "auto":
            template = (
                "{apa_b_executable} --bam-manifest {bam_manifest} --fasta {fasta} --gtf {gtf} "
                "--species {species} --assembly {assembly} --outdir {outdir} --threads {threads} "
                "--endpoint-source {endpoint_source} --endpoint-workers {endpoint_workers} "
                "--cluster-workers {cluster_workers} --deepip-threads {deepip_threads} "
                "--library-protocol {library_protocol} "
                "--validation-manifest {validation_manifest} --installation-manifest {installation_manifest}"
            )
        try:
            command = shlex.split(template.format(**replacements))
        except KeyError as exc:
            raise RuntimeError(f"Unknown placeholder in APA_B_COMMAND_TEMPLATE: {exc}") from exc
        run(command, log_dir / "apa_b" / genome / "engine.log", dry_run)
        if dry_run:
            continue
        catalog = genome_dir / "pas_catalog.tsv"
        counts = genome_dir / "pas_counts.tsv"
        deepip = genome_dir / "deepip_audit.tsv"
        engine_provenance = genome_dir / "engine_provenance.json"
        adapter_audit = genome_dir / "adapter_audit.json"
        _validate_polyaseqtrap_outputs(catalog, counts, deepip, [sample["sample_id"] for sample in samples])
        _validate_engine_provenance(engine_provenance, accepted_validation, genome)
        stats_dir = genome_dir / "drimseq"
        index = stats_dir / "result_index.tsv"
        genome_plan = RunPlan(plan.project, samples,
            [row for row in plan.sample_rows if row["genome"] == genome], contrasts, reference, {genome: reference})
        run_r_contrasts(
            module=f"apa_b_{genome}", plan=genome_plan, results=results,
            script=script_root / "R" / "drimseq_stager_all_pairs.R",
            common_arguments=["--counts", str(counts), "--catalog", str(catalog),
                "--samples", str(results / "00_metadata" / "validated_samples.tsv"),
                "--contrasts", str(results / "00_metadata" / "contrasts.tsv"),
                "--outdir", str(stats_dir), "--design", str(plan.project["design"]),
                "--fdr", str(plan.project["reporting"]["fdr"])],
            outdir=stats_dir, log_dir=log_dir / "apa_b" / genome / "contrasts",
            receipt_root=stats_dir / ".receipts", index_path=index,
            parallel_jobs=plan.project["resources"]["apa_b"]["contrast_parallel_jobs"],
            threads=1, memory_gb=plan.project["resources"]["apa_b"]["contrast_memory_gb"],
            output_suffixes=[
                ".drimseq_stager.tsv", ".gene_screen.tsv", ".gene_apa_summary.tsv",
                ".na_audit.tsv", ".fit_audit.tsv",
            ],
            signature_inputs=[counts, catalog], signature_parameters={
                "genome": genome, "design": plan.project["design"],
                "reporting": plan.project["reporting"], "settings": settings,
            },
            dry_run=False, force=force,
        )
        pcpa_catalog = genome_dir / "pcpa_candidate_catalog.tsv"
        pcpa_result = genome_dir / "candidate_pcpa.tsv"
        _classify_pcpa(catalog, reference["gtf"], pcpa_catalog)
        _filter_pcpa(pcpa_catalog, catalog, index, pcpa_result,
                     float(plan.project["reporting"]["fdr"]),
                     float(plan.project["reporting"]["min_abs_delta_pau"]))
        outputs.extend([
            manifest, catalog, counts, deepip, engine_provenance, adapter_audit, index,
            pcpa_catalog, pcpa_result,
        ])
        with index.open(encoding="utf-8", newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                _validate_na_audit(Path(row["na_audit_file"]))
                _validate_fit_audit(Path(row["fit_audit_file"]))
                outputs.extend([
                    Path(row["result_file"]), Path(row["gene_screen_file"]),
                    Path(row["gene_summary_file"]), Path(row["na_audit_file"]),
                    Path(row["fit_audit_file"]),
                ])
    if dry_run:
        event(log_dir, "apa_b", "dry_run", "Would run independent genome-specific APA-B and DRIMSeq/stageR")
        return
    write_receipt("apa_b", module_dir, signature, outputs, ["rna-ends2tracks", "apa-b"])
    event(log_dir, "apa_b", "completed", "Independent genome-specific APA-B + DRIMSeq/stageR branch")


def _classify_pcpa(catalog_path: Path, gtf: str, output: Path) -> None:
    genes = load_gene_models(gtf)
    with catalog_path.open(encoding="utf-8-sig", newline="") as handle:
        catalog = list(csv.DictReader(handle, delimiter="\t"))
    terminal_genes = {row["gene_id"] for row in catalog if row["gene_id"] and row["feature_class"].startswith("terminal")}
    rows: list[dict[str, str | int]] = []
    for row in catalog:
        gene_id = row["gene_id"]
        if gene_id not in terminal_genes or gene_id not in genes:
            continue
        feature = row["feature_class"]
        normalized_feature = feature.lower()
        if normalized_feature not in {"intron", "intronic", "internal_exon", "internal_exon_cds", "cds"}:
            continue
        position = int(row["start"]); gene = genes[gene_id]
        terminal = gene.end - 1 if gene.strand == "+" else gene.start
        if not (position < terminal if gene.strand == "+" else position > terminal):
            continue
        rows.append({
            "pas_id": row["pas_id"], "gene_id": gene_id, "chrom": row["chrom"], "start": position,
            "end": row["end"], "strand": row["strand"], "feature_class": feature,
            "consequence": "coding_truncating_intronic_PCPA" if normalized_feature in {"intron", "intronic", "internal_exon_cds", "cds"} else "upstream_exonic_termination",
            "interpretation": "candidate PCPA consistent with premature transcription termination",
        })
    headers = ["pas_id", "gene_id", "chrom", "start", "end", "strand", "feature_class", "consequence", "interpretation"]
    with output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=headers, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)


def _filter_pcpa(candidates_path: Path, catalog_path: Path, index_path: Path, output: Path, fdr: float, min_delta: float) -> None:
    with candidates_path.open(encoding="utf-8", newline="") as handle:
        candidates = {row["pas_id"]: row for row in csv.DictReader(handle, delimiter="\t")}
    with catalog_path.open(encoding="utf-8", newline="") as handle:
        catalog = {row["pas_id"]: row for row in csv.DictReader(handle, delimiter="\t")}
    selected: list[dict[str, str]] = []
    with index_path.open(encoding="utf-8", newline="") as handle:
        indexes = list(csv.DictReader(handle, delimiter="\t"))
    for index in indexes:
        with Path(index["result_file"]).open(encoding="utf-8", newline="") as handle:
            results = list(csv.DictReader(handle, delimiter="\t"))
            tested_terminal_genes = {catalog[row["feature_id"]]["gene_id"] for row in results
                                     if row.get("feature_id") in catalog and catalog[row["feature_id"]]["feature_class"].lower().startswith("terminal")}
            for row in results:
                adjusted = row.get("stageR_adjusted", "")
                if row.get("feature_id") not in candidates or adjusted in {"", "NA"}:
                    continue
                if candidates[row["feature_id"]]["gene_id"] not in tested_terminal_genes:
                    continue
                if float(adjusted) <= fdr and abs(float(row.get("delta_PAU", 0))) >= min_delta:
                    selected.append({**candidates[row["feature_id"]], "contrast_id": index["contrast_id"],
                                     "stageR_adjusted": adjusted, "delta_PAU": row.get("delta_PAU", "")})
    base = ["pas_id", "gene_id", "chrom", "start", "end", "strand", "feature_class", "consequence", "interpretation"]
    with output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=[*base, "contrast_id", "stageR_adjusted", "delta_PAU"], delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(selected)
