from __future__ import annotations

import csv
import gzip
import json
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

import pysam

from .config import RunPlan, sample_universe, signature_for
from .execution import run_bounded_processes
from .external import event, require_tools
from .mcell2019 import (
    assign_gene,
    build_gene_bins,
    count_sites,
    discover_sites,
    eligible_mapping,
    gene_counts,
    gtf_transcript_ends,
    internal_priming_at_position,
    load_gene_models,
    load_rescue_sites,
    pooled_cpm,
    read_chrom_sizes,
    rescue_overlap,
    transcript_end,
    write_tsv,
)
from .receipts import receipt_valid, write_receipt
from .statistics import run_r_contrasts

POSITION_FIELDS = ["chrom", "start", "end", "strand", "count"]


def _position_rows(counts: dict[tuple[str, str, int], int]) -> list[dict[str, int | str]]:
    return [
        {"chrom": chrom, "start": position, "end": position + 1, "strand": strand, "count": count}
        for (chrom, strand, position), count in sorted(counts.items())
    ]


def _write_positions(path: Path, counts: dict[tuple[str, str, int], int]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=POSITION_FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(_position_rows(counts))


def _read_positions(path: Path) -> dict[tuple[str, str, int], int]:
    result: dict[tuple[str, str, int], int] = {}
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            result[(row["chrom"], row["strand"], int(row["start"]))] = int(row["count"])
    return result


def _write_differential_pcpa(
    candidate_path: Path, index_path: Path, output: Path, fdr: float, min_delta: float,
) -> None:
    with candidate_path.open(encoding="utf-8", newline="") as handle:
        candidates = {row["pas_id"]: row for row in csv.DictReader(handle, delimiter="\t")}
    selected: list[dict[str, Any]] = []
    with index_path.open(encoding="utf-8", newline="") as handle:
        index_rows = list(csv.DictReader(handle, delimiter="\t"))
    for index in index_rows:
        with Path(index["result_file"]).open(encoding="utf-8", newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                pas_id = row.get("pas_id", "")
                adjusted = row.get("padj", "")
                delta = row.get("delta_PAU", "")
                if pas_id not in candidates or adjusted in {"", "NA"} or delta in {"", "NA"}:
                    continue
                if float(adjusted) <= fdr and abs(float(delta)) >= min_delta:
                    selected.append({**candidates[pas_id], "contrast_id": index["contrast_id"],
                                     "padj": adjusted, "delta_PAU": delta})
    candidate_fields = list(next(iter(candidates.values()))) if candidates else [
        "pas_id", "genome", "chrom", "summit_start", "summit_end", "start", "end",
        "actual_width", "strand", "gene_id", "feature_class", "assignment_status", "method",
        "interpretation",
    ]
    write_tsv(output, selected, [*candidate_fields, "contrast_id", "padj", "delta_PAU"])


def _sample_bams(plan: RunPlan, results: Path, sample_id: str) -> tuple[list[Path], str]:
    lanes = [
        results / "02_alignment" / row["sample_id"] / "lanes" /
        f"{row['sample_id']}.{row['technical_replicate_id']}.{row['lane_id']}.bam"
        for row in plan.sample_rows if row["sample_id"] == sample_id
    ]
    present = [path for path in lanes if path.is_file()]
    if len(present) == len(lanes) and present:
        return present, "all_alignment_records"
    final = results / "02_alignment" / sample_id / f"{sample_id}.bam"
    if final.is_file():
        return [final], "C0_only_after_intermediate_cleanup"
    raise RuntimeError(f"No alignment BAMs available for exact-end extraction: {sample_id}")


def _extract_sample(plan: RunPlan, results: Path, sample: dict[str, str], force: bool) -> tuple[Path, Path, Path, Path, Path]:
    sample_id, genome = sample["sample_id"], sample["genome"]
    reference = plan.reference_for(genome)
    outdir = results / "03_exact_ends" / genome / sample_id
    receipt_dir = outdir / ".receipt"
    c1_path = outdir / "C1_exact_ends.tsv.gz"
    c1s_path = outdir / "C1S_uncertain_ends.tsv.gz"
    c2_path = outdir / "C2_filtered_ends.tsv.gz"
    c2r_path = outdir / "C2R_internal_priming_rejects.tsv.gz"
    audit_path = outdir / "end_audit.json"
    bams, audit_scope = _sample_bams(plan, results, sample_id)
    settings = plan.project["apa_a"]
    signature = signature_for([*bams, reference["fasta"], reference["gtf"]], {
        "module": "exact_ends", "sample": sample, "settings": settings,
        "atlas": reference.get("pas_atlas", ""),
    })
    outputs = [c1_path, c1s_path, c2_path, c2r_path, audit_path]
    if not force and receipt_valid(receipt_dir, signature):
        return tuple(outputs)  # type: ignore[return-value]

    c1: dict[tuple[str, str, int], int] = defaultdict(int)
    c1s: dict[tuple[str, str, int], int] = defaultdict(int)
    audit: Counter[str] = Counter()
    compatibility = settings.get("mapping_policy") == "legacy_random_multimapper"
    for bam_path in bams:
        with pysam.AlignmentFile(bam_path, "rb") as bam:
            for read in bam.fetch(until_eof=True):
                eligible, mapping_class = eligible_mapping(read, compatibility)
                audit[mapping_class] += 1
                if read.is_duplicate:
                    audit["duplicate_flagged"] += 1
                if not eligible:
                    continue
                position, strand, clipped = transcript_end(read)
                key = (read.reference_name, strand, position)
                if clipped:
                    c1s[key] += 1
                    audit["end_defining_clipped"] += 1
                else:
                    c1[key] += 1

    fasta = pysam.FastaFile(reference["fasta"])
    rescue_atlas = load_rescue_sites(reference.get("pas_atlas"), settings["mask_rescue_tier"])
    annotated_ends = gtf_transcript_ends(reference["gtf"])
    c2: dict[tuple[str, str, int], int] = {}
    c2r: dict[tuple[str, str, int], int] = {}
    for key, value in c1.items():
        chrom, strand, position = key
        masked = internal_priming_at_position(
            fasta, chrom, position, strand,
            int(settings["internal_priming_consecutive_bases"]),
            int(settings["internal_priming_window_nt"]),
            int(settings["internal_priming_min_bases_in_window"]),
        )
        rescued = masked and rescue_overlap(chrom, strand, position, (annotated_ends, rescue_atlas), 20)
        (c2 if not masked or rescued else c2r)[key] = value
        audit["mask_rescued" if rescued else "mask_rejected" if masked else "mask_passed"] += value
    fasta.close()
    audit.update({
        "C0": sum(c1.values()) + sum(c1s.values()), "C1": sum(c1.values()), "C1S": sum(c1s.values()),
        "C2": sum(c2.values()), "C2R": sum(c2r.values()), "duplicate_policy_retained": 1,
    })
    if audit["C0"] != audit["C1"] + audit["C1S"] or audit["C1"] != audit["C2"] + audit["C2R"]:
        raise RuntimeError(f"C0-C2R count invariant failed for {sample_id}")
    outdir.mkdir(parents=True, exist_ok=True)
    for path, counts in ((c1_path, c1), (c1s_path, c1s), (c2_path, c2), (c2r_path, c2r)):
        _write_positions(path, counts)
    audit_payload = {**dict(sorted(audit.items())), "sample_id": sample_id, "genome": genome, "audit_scope": audit_scope}
    audit_path.write_text(json.dumps(audit_payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_receipt("exact_ends_sample", receipt_dir, signature, outputs, ["rna-ends2tracks", "exact_ends", sample_id])
    return tuple(outputs)  # type: ignore[return-value]


def _extract_sample_process(
    plan: RunPlan, results: Path, sample: dict[str, str], force: bool,
) -> tuple[Path, Path, Path, Path, Path]:
    """Pickle-safe process-pool entry point for one biological sample."""
    return _extract_sample(plan, results, sample, force)


def _discover_genome(plan: RunPlan, results: Path, genome: str) -> tuple[list[Path], list[dict[str, Any]]]:
    reference = plan.reference_for(genome)
    samples = [sample for sample in plan.samples if sample["genome"] == genome]
    sample_ids = [sample["sample_id"] for sample in samples]
    sample_counts = {
        sample_id: _read_positions(results / "03_exact_ends" / genome / sample_id / "C2_filtered_ends.tsv.gz")
        for sample_id in sample_ids
    }
    pooled, totals = pooled_cpm(sample_counts)
    output = results / "04_active_pas" / genome
    output.mkdir(parents=True, exist_ok=True)
    pooled_path = output / "condition_blind_pooled_C2_CPM.tsv.gz"
    with gzip.open(pooled_path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["chrom", "start", "end", "strand", "pooled_cpm"])
        for (chrom, strand, position), value in sorted(pooled.items()):
            writer.writerow([chrom, position, position + 1, strand, f"{value:.12g}"])
    settings = plan.project["apa_a"]
    sites = discover_sites(
        pooled, read_chrom_sizes(reference["chrom_sizes"]),
        int(settings["discovery_window_nt"]), float(settings["discovery_threshold"]),
        int(settings["discovery_rounds"]),
    )
    if not sites:
        raise RuntimeError(
            f"No active PAS passed Mcell2019 discovery for {genome}; inspect C2 totals and PAS thresholds"
        )
    genes = load_gene_models(reference["gtf"])
    gene_bins = build_gene_bins(genes, int(settings["gene_downstream_extension_nt"]))
    catalog: list[dict[str, Any]] = []
    for index, site in enumerate(sites, start=1):
        gene_id, feature, assignment = assign_gene(
            site, genes, int(settings["gene_downstream_extension_nt"]), gene_bins
        )
        catalog.append({
            "pas_id": f"{genome}_PAS_{index:08d}", "genome": genome, "chrom": site.chrom,
            "summit_start": site.summit, "summit_end": site.summit + 1,
            "start": site.start, "end": site.end, "actual_width": site.end - site.start,
            "strand": site.strand, "gene_id": gene_id, "feature_class": feature,
            "assignment_status": assignment, "method": "mcell2019_two_round_30nt",
        })
    counts = count_sites(sites, sample_counts)
    for row, site_row in zip(counts, catalog):
        row["pas_id"] = site_row["pas_id"]
    c4 = gene_counts(catalog, counts, sample_ids)
    catalog_path, counts_path, c4_path = output / "active_pas_catalog.tsv", output / "C3_active_pas_counts.tsv", output / "C4_active_pas_gene_counts.tsv"
    write_tsv(catalog_path, catalog)
    write_tsv(counts_path, counts, ["pas_id", *sample_ids])
    write_tsv(c4_path, c4, ["gene_id", *sample_ids])
    write_tsv(output / "ambiguous_pas.tsv", [row for row in catalog if row["assignment_status"] == "ambiguous_multi_gene"], list(catalog[0]) if catalog else [])
    write_tsv(output / "C2_library_totals.tsv", [{"sample_id": key, "C2_total": value} for key, value in totals.items()])
    terminal_genes = {row["gene_id"] for row in catalog
                      if row["assignment_status"] == "unique" and row["feature_class"] == "terminal_exon"}
    pcpa_candidates = [
        {**row, "interpretation": "candidate intragenic premature cleavage/polyadenylation site"}
        for row in catalog
        if row["assignment_status"] == "unique" and row["gene_id"] in terminal_genes
        and row["feature_class"] in {"intron", "other_exon"}
    ]
    pcpa_path = output / "pcpa_candidate_catalog.tsv"
    write_tsv(pcpa_path, pcpa_candidates, [*list(catalog[0]), "interpretation"])
    c2_sum = sum(sum(values.values()) for values in sample_counts.values())
    c3_sum = sum(sum(int(row[sample]) for sample in sample_ids) for row in counts)
    if c3_sum > c2_sum:
        raise RuntimeError(f"C3 count invariant failed for {genome}: {c3_sum} > C2 {c2_sum}")
    audit_path = output / "count_universe_audit.json"
    audit_path.write_text(json.dumps({
        "genome": genome, "samples": len(sample_ids), "active_pas": len(catalog),
        "C2_total": c2_sum, "C3_total": c3_sum, "C3_le_C2": True,
        "nonoverlapping_pas_intervals": True, "one_C2_end_per_C3_maximum": True,
        "ambiguous_pas": sum(row["assignment_status"] == "ambiguous_multi_gene" for row in catalog),
        "pcpa_candidate_sites": len(pcpa_candidates),
    }, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return [pooled_path, catalog_path, counts_path, c4_path, output / "ambiguous_pas.tsv",
            output / "C2_library_totals.tsv", pcpa_path, audit_path], catalog


def exact_ends_stage(plan: RunPlan, results: Path, dry_run: bool = False, force: bool = False) -> None:
    log_dir = results / "logs"
    if dry_run:
        event(log_dir, "exact_ends", "dry_run", "Would build C1/C1S/C2/C2R for every biological sample")
        return
    exact_root = results / "03_exact_ends"
    signature_inputs = [
        results / "02_alignment" / sample["sample_id"] / f"{sample['sample_id']}.bam" for sample in plan.samples
    ]
    exact_signature = signature_for(signature_inputs, {
        "module": "exact_ends", "samples": plan.samples, "apa_a": plan.project["apa_a"],
        "references": plan.references or {plan.reference["assembly"]: plan.reference},
    })
    if not force and receipt_valid(exact_root, exact_signature):
        event(log_dir, "exact_ends", "skipped", "Valid project receipt")
    else:
        jobs = [
            (sample["sample_id"], _extract_sample_process, (plan, results, sample, force))
            for sample in plan.samples
        ]
        parallel_jobs = plan.project["resources"]["apa_a"]["extraction_parallel_jobs"]
        event(log_dir, "exact_ends", "started",
              f"Launching {len(jobs)} sample workers with process_parallel_jobs={parallel_jobs}")
        outputs_nested = run_bounded_processes(
            "exact_ends", jobs, plan.project["resources"]["apa_a"]["extraction_parallel_jobs"],
            results / ".checkpoints" / "timings" / "exact_ends",
            progress=lambda label, status: event(
                log_dir, "exact_ends", status, f"Sample worker {label} {status}"
            ),
        )
        outputs = [path for paths in outputs_nested for path in paths]
        write_receipt("exact_ends", exact_root, exact_signature, outputs, ["rna-ends2tracks", "exact_ends"])
    event(log_dir, "exact_ends", "completed", f"Validated C0-C2R invariants for {len(plan.samples)} samples")


def active_pas_stage(plan: RunPlan, results: Path, dry_run: bool = False, force: bool = False) -> None:
    log_dir = results / "logs"
    if dry_run:
        event(log_dir, "active_pas", "dry_run", "Would pool C2 CPM condition-blind and run two Mcell2019 discovery rounds")
        return
    active_root = results / "04_active_pas"
    sample_set = sample_universe(plan)
    universe_path = active_root / "sample_set.json"
    if universe_path.is_file():
        prior = json.loads(universe_path.read_text(encoding="utf-8"))
        if prior != sample_set:
            raise RuntimeError(
                "The active-PAS sample universe changed. Use a new OUTPUT_DIR; condition-blind PAS catalogs "
                "must never be silently reused or overwritten after sample membership changes."
            )
    c2_paths = [
        results / "03_exact_ends" / sample["genome"] / sample["sample_id"] / "C2_filtered_ends.tsv.gz"
        for sample in plan.samples
    ]
    active_signature = signature_for(c2_paths, {
        "module": "active_pas", "sample_set": sample_set, "settings": plan.project["apa_a"],
    })
    if not force and receipt_valid(active_root, active_signature):
        event(log_dir, "active_pas", "skipped", "Valid project-wide PAS-universe receipt")
        return

    active_outputs: list[Path] = []
    for genome in plan.references or {plan.reference["assembly"]: plan.reference}:
        if any(sample["genome"] == genome for sample in plan.samples):
            outputs, _ = _discover_genome(plan, results, genome)
            active_outputs.extend(outputs)
    active_root.mkdir(parents=True, exist_ok=True)
    universe_path.write_text(json.dumps(sample_set, indent=2) + "\n", encoding="utf-8")
    active_outputs.append(universe_path)
    write_receipt("active_pas", active_root, active_signature, active_outputs, ["rna-ends2tracks", "active_pas"])
    event(log_dir, "active_pas", "completed", "Created genome-specific condition-blind active-PAS universes")


def apa_statistics_stage(
    plan: RunPlan, results: Path, script_root: Path, dry_run: bool = False, force: bool = False,
) -> None:
    log_dir = results / "logs"
    if dry_run:
        event(log_dir, "apa_a_statistics", "dry_run", "Would run per-contrast DEXSeq and Mcell2019 shift classification")
        return
    require_tools(["Rscript"])
    active_root = results / "04_active_pas"
    stats_root = results / "06_apa_a_mcell2019"
    stats_outputs: list[Path] = []
    active_inputs: list[Path] = []
    for genome in plan.references or {plan.reference["assembly"]: plan.reference}:
        genome_samples = [sample for sample in plan.samples if sample["genome"] == genome]
        genome_contrasts = [contrast for contrast in plan.contrasts if contrast["genome"] == genome]
        if not genome_contrasts:
            continue
        active_inputs.extend([
            active_root / genome / "C3_active_pas_counts.tsv",
            active_root / genome / "active_pas_catalog.tsv",
        ])
        genome_plan = RunPlan(plan.project, genome_samples, [row for row in plan.sample_rows if row["genome"] == genome], genome_contrasts, plan.reference_for(genome), {genome: plan.reference_for(genome)})
        outdir = stats_root / genome / "dexseq"
        index_path = outdir / "result_index.tsv"
        run_r_contrasts(
            module=f"apa_a_{genome}", plan=genome_plan, results=results,
            script=script_root / "R" / "dexseq_all_pairs.R",
            common_arguments=[
                "--counts", str(active_root / genome / "C3_active_pas_counts.tsv"),
                "--catalog", str(active_root / genome / "active_pas_catalog.tsv"),
                "--samples", str(results / "00_metadata" / "validated_samples.tsv"),
                "--contrasts", str(results / "00_metadata" / "contrasts.tsv"),
                "--outdir", str(outdir), "--min-count", "5", "--design", str(plan.project["design"]),
                "--fdr", str(plan.project["reporting"]["fdr"]),
            ],
            outdir=outdir, log_dir=results / "logs" / "apa_a" / genome,
            receipt_root=outdir / ".receipts", index_path=index_path,
            parallel_jobs=plan.project["resources"]["apa_a"]["contrast_parallel_jobs"], threads=1,
            memory_gb=plan.project["resources"]["apa_a"]["contrast_memory_gb"],
            output_suffixes=[".dexseq.tsv", ".apa_shift.tsv"],
            signature_inputs=[active_root / genome / "C3_active_pas_counts.tsv", active_root / genome / "active_pas_catalog.tsv"],
            signature_parameters={
                "design": plan.project["design"], "genome": genome,
                "reporting": plan.project["reporting"],
            }, dry_run=False, force=force,
        )
        pcpa_path = stats_root / genome / "candidate_pcpa.tsv"
        _write_differential_pcpa(
            active_root / genome / "pcpa_candidate_catalog.tsv", index_path, pcpa_path,
            float(plan.project["reporting"]["fdr"]),
            float(plan.project["reporting"]["min_abs_delta_pau"]),
        )
        stats_outputs.extend([index_path, pcpa_path])
        with index_path.open(encoding="utf-8", newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                stats_outputs.extend([Path(row["result_file"]), Path(row["shift_file"])])
    if stats_outputs:
        stats_signature = signature_for(active_inputs, {
            "module": "apa_a_statistics", "contrasts": plan.contrasts,
            "design": plan.project["design"], "reporting": plan.project["reporting"],
        })
        write_receipt("apa_a_mcell2019", stats_root, stats_signature, stats_outputs, ["rna-ends2tracks", "apa-a"])
    event(log_dir, "apa_a_statistics", "completed", "Completed genome-specific DEXSeq and shift analyses")


def apa_a_mcell2019(
    plan: RunPlan, results: Path, script_root: Path, dry_run: bool = False, force: bool = False,
) -> None:
    exact_ends_stage(plan, results, dry_run, force)
    active_pas_stage(plan, results, dry_run, force)
    apa_statistics_stage(plan, results, script_root, dry_run, force)
