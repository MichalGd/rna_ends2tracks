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


def _validation_manifest(path: Path, genomes: list[str]) -> dict[str, object]:
    if not path.is_file():
        raise RuntimeError(f"APA-B validation manifest is unavailable: {path}")
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError) as exc:
        raise RuntimeError(f"APA-B validation manifest is invalid JSON: {path}") from exc
    if payload.get("schema_version") != 1 or payload.get("status") != "accepted":
        raise RuntimeError("APA-B validation manifest must have schema_version=1 and status=accepted")
    engine = payload.get("engine", {})
    model = payload.get("model", {})
    models = payload.get("models", {})
    environment = payload.get("environment", {})
    if not isinstance(engine, dict) or not re.fullmatch(r"[0-9a-fA-F]{7,40}", str(engine.get("source_commit", ""))):
        raise RuntimeError("APA-B validation manifest requires a pinned engine.source_commit")
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
    if "quantseq_rev_v2_se" not in payload.get("library_protocols", []):
        raise RuntimeError("APA-B validation does not cover quantseq_rev_v2_se")
    missing = sorted(set(genomes).difference(map(str, payload.get("assemblies", []))))
    if missing:
        raise RuntimeError("APA-B validation does not cover assemblies: " + ", ".join(missing))
    pilot = payload.get("pilot", {})
    if not isinstance(pilot, dict) or pilot.get("synthetic_pass") is not True:
        raise RuntimeError("APA-B validation requires a passing synthetic coordinate/strand pilot")
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
    for section, field in (("engine", "source_commit"), ("model", "sha256"), ("environment", "sha256")):
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
    accepted_validation = _validation_manifest(validation_path, list(plan.references))
    bams = [results / "02_alignment" / sample["sample_id"] / f"{sample['sample_id']}.bam" for sample in plan.samples]
    ref_inputs = [Path(plan.reference_for(genome)[key]) for genome in plan.references for key in ("fasta", "gtf")]
    signature = signature_for([*bams, *ref_inputs, validation_path, installation_path], {
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
            writer.writerow(["sample_id", "bam"])
            writer.writerows((sample["sample_id"], bam) for sample, bam in zip(samples, genome_bams))
        replacements = {"bam_manifest": str(manifest), "fasta": reference["fasta"],
                        "gtf": reference["gtf"], "outdir": str(genome_dir),
                        "species": reference["species"], "assembly": reference["assembly"],
                        "threads": str(plan.project["resources"]["apa_b"]["engine_threads"]),
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
            output_suffixes=[".drimseq_stager.tsv", ".gene_screen.tsv", ".gene_apa_summary.tsv"],
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
                outputs.extend([
                    Path(row["result_file"]), Path(row["gene_screen_file"]),
                    Path(row["gene_summary_file"]),
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
