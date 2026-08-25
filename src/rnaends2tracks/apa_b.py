from __future__ import annotations

import csv
import shlex
from pathlib import Path

from .apa_a import load_genes
from .config import RunPlan, signature_for
from .external import event, require_tools, run
from .receipts import receipt_valid, write_receipt


REQUIRED_CATALOG = {"pas_id", "gene_id", "chrom", "start", "end", "strand", "feature_class"}


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
    module_dir = results / "05_apa_b_polyaseqtrap_drimseq"
    log_dir = results / "provenance" / "logs"
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
    bams = [results / "02_alignment" / sample["sample_id"] / f"{sample['sample_id']}.bam" for sample in plan.samples]
    if not dry_run and any(not path.is_file() for path in bams):
        raise RuntimeError("APA-B requires all shared BAMs")
    signature = signature_for([*bams, plan.reference["fasta"], plan.reference["gtf"]], {
        "module": "apa_b", "settings": settings, "design": plan.project["design"],
        "samples": plan.samples, "contrasts": plan.contrasts, "reporting": plan.project.get("reporting", {}),
    }) if not dry_run else "dry-run"
    catalog = module_dir / "pas_catalog.tsv"
    counts = module_dir / "pas_counts.tsv"
    deepip = module_dir / "deepip_audit.tsv"
    stats_index = module_dir / "drimseq" / "result_index.tsv"
    pcpa_catalog = module_dir / "pcpa_candidate_catalog.tsv"
    pcpa_result = module_dir / "candidate_pcpa.tsv"
    if not force and not dry_run and receipt_valid(module_dir, signature):
        event(log_dir, "apa_b", "skipped", "Valid matching receipt")
        return
    manifest = module_dir / "bam_manifest.tsv"
    with manifest.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["sample_id", "bam"])
        writer.writerows((sample["sample_id"], bam) for sample, bam in zip(plan.samples, bams))
    replacements = {
        "bam_manifest": str(manifest), "fasta": plan.reference["fasta"], "gtf": plan.reference["gtf"],
        "outdir": str(module_dir), "species": plan.reference["species"], "assembly": plan.reference["assembly"],
    }
    try:
        command = shlex.split(template.format(**replacements))
    except KeyError as exc:
        raise RuntimeError(f"Unknown placeholder in apa_b.command_template: {exc}") from exc
    run(command, log_dir / "apa_b" / "polyaseqtrap.log", dry_run)
    if dry_run:
        event(log_dir, "apa_b", "dry_run", "Would run independent PolyAseqTrap adapter and DRIMSeq/stageR")
        return
    _validate_polyaseqtrap_outputs(catalog, counts, deepip, [sample["sample_id"] for sample in plan.samples])
    require_tools(["Rscript"])
    run([
        "Rscript", str(script_root / "R" / "drimseq_stager_all_pairs.R"),
        "--counts", str(counts), "--catalog", str(catalog),
        "--samples", str(results / "00_metadata" / "validated_samples.tsv"),
        "--contrasts", str(results / "00_metadata" / "contrasts.tsv"),
        "--outdir", str(module_dir / "drimseq"),
        "--design", str(plan.project["design"]),
        "--fdr", str(plan.project.get("reporting", {}).get("fdr", 0.05)),
    ], log_dir / "apa_b" / "drimseq_stager.log", False)
    _classify_pcpa(catalog, plan.reference["gtf"], pcpa_catalog)
    _filter_pcpa(
        pcpa_catalog, catalog, stats_index, pcpa_result,
        float(plan.project.get("reporting", {}).get("fdr", 0.05)),
        float(plan.project.get("reporting", {}).get("min_abs_delta_pau", 0.10)),
    )
    write_receipt("apa_b", module_dir, signature, [catalog, counts, deepip, stats_index, pcpa_catalog, pcpa_result], ["rna-ends2tracks", "apa-b"])
    event(log_dir, "apa_b", "completed", "Independent PolyAseqTrap + DRIMSeq/stageR branch")


def _classify_pcpa(catalog_path: Path, gtf: str, output: Path) -> None:
    genes, _ = load_genes(gtf)
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
