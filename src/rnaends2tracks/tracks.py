from __future__ import annotations

import csv
import gzip
import json
from pathlib import Path
from typing import Any

from .config import RunPlan, signature_for
from .execution import run_bounded_processes
from .external import event, progress_events, require_tools, run, run_capture, run_to_path
from .paths import workflow_asset
from .receipts import receipt_valid, write_receipt

END_FAMILIES = {
    "exact_ends": "C1_exact_ends.tsv.gz",
    "filtered_ends": "C2_filtered_ends.tsv.gz",
    "rejected_ends": "C2R_internal_priming_rejects.tsv.gz",
}
TRACK_FIELDS = ["sample_id", "genome", "family", "normalization", "denominator", "scale", "count_universe"]


def _negated(source: Path, target: Path) -> None:
    with source.open(encoding="utf-8") as inp, target.open("w", encoding="utf-8") as out:
        for line in inp:
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 4:
                fields[3] = str(-float(fields[3]))
            out.write("\t".join(fields) + "\n")


def _end_counts(path: Path) -> list[dict[str, str]]:
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _active_counts(catalog: Path, matrix: Path, sample_id: str) -> tuple[list[dict[str, str]], int]:
    with catalog.open(encoding="utf-8", newline="") as handle:
        metadata = {row["pas_id"]: row for row in csv.DictReader(handle, delimiter="\t")}
    rows: list[dict[str, str]] = []
    assigned_total = 0
    with matrix.open(encoding="utf-8", newline="") as handle:
        for count_row in csv.DictReader(handle, delimiter="\t"):
            meta = metadata[count_row["pas_id"]]
            count = int(count_row[sample_id])
            if meta["assignment_status"] == "unique":
                assigned_total += count
            rows.append({
                "chrom": meta["chrom"], "start": meta["summit_start"],
                "end": str(int(meta["summit_start"]) + 1), "strand": meta["strand"], "count": str(count),
            })
    return rows, assigned_total


def _write_signal_bedgraphs(
    rows: list[dict[str, str]], plus: Path, minus: Path, scale: float,
    chromosome_order: dict[str, int],
) -> None:
    plus.parent.mkdir(parents=True, exist_ok=True)
    ordered = sorted(rows, key=lambda row: (
        chromosome_order.get(row["chrom"], len(chromosome_order)), int(row["start"]), int(row["end"])
    ))
    with plus.open("w", encoding="utf-8") as plus_handle, minus.open("w", encoding="utf-8") as minus_handle:
        for row in ordered:
            value = int(row["count"]) * scale
            handle = plus_handle if row["strand"] == "+" else minus_handle
            display = value if row["strand"] == "+" else -value
            handle.write(f"{row['chrom']}\t{row['start']}\t{row['end']}\t{display:.12g}\n")


def _size_factors(path: Path) -> dict[str, dict[str, float]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return {
            row["sample_id"]: {
                "deseq2": float(row["deseq2_scale"]), "robust_cpm": float(row["robust_cpm_scale"]),
                "size_factor": float(row["size_factor"]),
            }
            for row in csv.DictReader(handle, delimiter="\t")
        }


def _bam_count(bam: Path, log: Path, threads: int) -> int:
    stdout = run_capture(["samtools", "view", "-@", str(threads), "-c", str(bam)], log)
    return int(stdout.strip())


def _strand_bam(bam: Path, strand: str, output: Path, threads: int, log: Path) -> None:
    flag_args = ["-f", "16"] if strand == "plus" else ["-F", "16"]
    run(["samtools", "view", "-@", str(threads), "-b", *flag_args, "-o", str(output), str(bam)], log, False)


def _all_read_bedgraph(strand_bam: Path, strand: str, output: Path, scale: float, log: Path) -> None:
    positive = output if strand == "plus" else output.with_suffix(".positive.bedGraph")
    run_to_path(
        ["bedtools", "genomecov", "-ibam", str(strand_bam), "-bg", "-scale", f"{scale:.15g}"],
        positive, log,
    )
    if strand == "minus":
        _negated(positive, output)


def _sample_tracks_subset(
    plan: RunPlan,
    results: Path,
    sample: dict[str, str],
    force: bool,
    selected_families: tuple[str, ...],
    selected_normalizations: tuple[str, ...],
    receipt_group: str,
) -> tuple[list[Path], list[dict[str, Any]]]:
    sample_id, genome = sample["sample_id"], sample["genome"]
    reference = plan.reference_for(genome)
    settings = plan.project["tracks"]
    families = {
        key: bool(value and key in selected_families)
        for key, value in settings["families"].items()
    }
    enabled_normalizations = {
        key: bool(value and key in selected_normalizations)
        for key, value in settings["normalizations"].items()
    }
    resource = plan.project["resources"]["tracks"]
    track_root = results / "09_tracks"
    intermediate = track_root / ".intermediate" / sample_id
    receipt_dir = track_root / ".receipts" / receipt_group / sample_id
    bam = results / "02_alignment" / sample_id / f"{sample_id}.bam"
    size_factor_path = results / "04_active_pas" / genome / "C4_track_size_factors.tsv"
    inputs: list[Path] = [bam, Path(reference["chrom_sizes"])]
    needs_factors = (
        enabled_normalizations["deseq2"] or enabled_normalizations["robust_cpm"]
    ) and (families["active_pas"] or families["filtered_ends"])
    if needs_factors:
        inputs.append(size_factor_path)
    for family, filename in END_FAMILIES.items():
        if families[family]:
            inputs.append(results / "03_exact_ends" / genome / sample_id / filename)
    if families["active_pas"]:
        inputs.extend([
            results / "04_active_pas" / genome / "active_pas_catalog.tsv",
            results / "04_active_pas" / genome / "C3_active_pas_counts.tsv",
        ])
    subset_settings = {
        "families": families,
        "normalizations": enabled_normalizations,
        "generate_bigwigs": settings["generate_bigwigs"],
        "retain_bedgraph": settings["retain_bedgraph"],
    }
    signature = signature_for(
        inputs,
        {"module": receipt_group, "sample": sample, "settings": subset_settings},
    )
    if not force and receipt_valid(receipt_dir, signature):
        receipt = receipt_dir / "run_receipt.json"
        fragment = receipt_dir / "normalization.tsv"
        receipt_payload = json.loads(receipt.read_text(encoding="utf-8"))
        with fragment.open(encoding="utf-8", newline="") as handle:
            return [Path(row["path"]) for row in receipt_payload["outputs"]], list(csv.DictReader(handle, delimiter="\t"))
    factors = _size_factors(size_factor_path)[sample_id] if needs_factors else {}
    expected: list[Path] = []
    normalization_rows: list[dict[str, Any]] = []
    intermediate.mkdir(parents=True, exist_ok=True)
    with Path(reference["chrom_sizes"]).open(encoding="utf-8") as handle:
        chromosome_sizes = [line.split()[:2] for line in handle if line.strip()]
    chromosome_order = {fields[0]: index for index, fields in enumerate(chromosome_sizes)}
    if not chromosome_sizes:
        raise RuntimeError(f"Chromosome-size file is empty: {reference['chrom_sizes']}")

    family_rows: dict[str, list[dict[str, str]]] = {}
    denominators: dict[str, int] = {}
    if families["all_reads"]:
        denominators["all_reads"] = _bam_count(bam, results / "logs" / "tracks" / f"{sample_id}.count.log", resource["samtools_threads"])
    for family, filename in END_FAMILIES.items():
        if families[family]:
            rows = _end_counts(results / "03_exact_ends" / genome / sample_id / filename)
            family_rows[family] = rows
            denominators[family] = sum(int(row["count"]) for row in rows)
    if families["rejected_ends"]:
        if "exact_ends" not in denominators:
            c1_rows = _end_counts(results / "03_exact_ends" / genome / sample_id / END_FAMILIES["exact_ends"])
            denominators["exact_ends"] = sum(int(row["count"]) for row in c1_rows)
        denominators["rejected_ends"] = denominators["exact_ends"]
    if families["active_pas"]:
        rows, assigned = _active_counts(
            results / "04_active_pas" / genome / "active_pas_catalog.tsv",
            results / "04_active_pas" / genome / "C3_active_pas_counts.tsv", sample_id,
        )
        family_rows["active_pas"] = rows
        denominators["active_pas"] = assigned

    strand_bams: dict[str, Path] = {}
    if families["all_reads"]:
        for strand in ("plus", "minus"):
            strand_bam = intermediate / f"{sample_id}.all_reads.{strand}.strand.bam"
            _strand_bam(
                bam, strand, strand_bam, resource["samtools_threads"],
                results / "logs" / "tracks" / f"{sample_id}.all_reads.{strand}.strand_bam.log",
            )
            strand_bams[strand] = strand_bam

    for family, enabled in families.items():
        if not enabled:
            continue
        normalization_scales: list[tuple[str, float]] = []
        if enabled_normalizations["raw"]:
            normalization_scales.append(("raw", 1.0))
        if enabled_normalizations["cpm"]:
            denominator = denominators.get(family, 0)
            normalization_scales.append(("cpm", 1_000_000.0 / denominator if denominator else 0.0))
        if family in {"filtered_ends", "active_pas"}:
            if enabled_normalizations["deseq2"]:
                normalization_scales.append(("deseq2", factors["deseq2"]))
            if enabled_normalizations["robust_cpm"]:
                normalization_scales.append(("robust_cpm", factors["robust_cpm"]))
        for normalization, scale in normalization_scales:
            normalization_rows.append({
                "sample_id": sample_id, "genome": genome, "family": family,
                "normalization": normalization, "denominator": denominators.get(family, ""),
                "scale": f"{scale:.15g}", "count_universe": {
                    "all_reads": "C0", "exact_ends": "C1", "filtered_ends": "C2",
                    "rejected_ends": "C1", "active_pas": "assigned_C3",
                }[family],
            })
            family_dir = track_root / family / normalization
            for strand in ("plus", "minus"):
                bedgraph = intermediate / f"{sample_id}.{family}.{normalization}.{strand}.bedGraph"
                bigwig = family_dir / f"{sample_id}.{family}.{normalization}.transcript_{strand}.bw"
                bigwig.parent.mkdir(parents=True, exist_ok=True)
                if family == "all_reads":
                    _all_read_bedgraph(
                        strand_bams[strand], strand, bedgraph, scale,
                        results / "logs" / "tracks" / f"{sample_id}.{family}.{normalization}.{strand}.log",
                    )
                else:
                    plus = intermediate / f"{sample_id}.{family}.{normalization}.plus.bedGraph"
                    minus = intermediate / f"{sample_id}.{family}.{normalization}.minus.bedGraph"
                    _write_signal_bedgraphs(family_rows[family], plus, minus, scale, chromosome_order)
                    bedgraph = plus if strand == "plus" else minus
                if settings["generate_bigwigs"]:
                    temporary_bigwig = bigwig.with_name(f".{bigwig.name}.tmp")
                    conversion_source = bedgraph
                    if bedgraph.stat().st_size == 0:
                        # UCSC bedGraphToBigWig rejects an empty input. A single
                        # zero-valued interval produces an honest zero-signal
                        # BigWig while the retained bedGraph, if requested,
                        # remains empty.
                        conversion_source = bedgraph.with_name(f".{bedgraph.name}.zero.bedGraph")
                        conversion_source.write_text(
                            f"{chromosome_sizes[0][0]}\t0\t{chromosome_sizes[0][1]}\t0\n",
                            encoding="utf-8",
                        )
                    run(["bedGraphToBigWig", str(conversion_source), reference["chrom_sizes"], str(temporary_bigwig)],
                        results / "logs" / "tracks" / f"{sample_id}.{family}.{normalization}.{strand}.log", False)
                    temporary_bigwig.replace(bigwig)
                    if conversion_source != bedgraph:
                        conversion_source.unlink()
                    expected.append(bigwig)
                if settings["retain_bedgraph"]:
                    retained = family_dir / bedgraph.name
                    temporary_retained = retained.with_name(f".{retained.name}.tmp")
                    temporary_retained.write_bytes(bedgraph.read_bytes())
                    temporary_retained.replace(retained)
                    expected.append(retained)
    receipt_dir.mkdir(parents=True, exist_ok=True)
    fragment = receipt_dir / "normalization.tsv"
    with fragment.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=TRACK_FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(normalization_rows)
    expected.append(fragment)
    write_receipt(
        receipt_group + "_sample", receipt_dir, signature, expected,
        ["rna-ends2tracks", receipt_group, sample_id],
    )
    return expected, normalization_rows


def _write_normalization(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp")
    with temporary.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=TRACK_FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


def _read_normalization(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        return []
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _receipt_outputs(receipt_dir: Path) -> list[Path]:
    receipt = receipt_dir / "run_receipt.json"
    if not receipt.is_file():
        return []
    payload = json.loads(receipt.read_text(encoding="utf-8"))
    return [Path(row["path"]) for row in payload.get("outputs", [])]


def _run_track_subset(
    plan: RunPlan,
    results: Path,
    force: bool,
    families: tuple[str, ...],
    normalizations: tuple[str, ...],
    receipt_group: str,
    event_module: str,
) -> tuple[list[Path], list[dict[str, Any]]]:
    jobs = [
        (
            sample["sample_id"],
            _sample_tracks_subset,
            (plan, results, sample, force, families, normalizations, receipt_group),
        )
        for sample in plan.samples
    ]
    completed = run_bounded_processes(
        receipt_group, jobs, plan.project["resources"]["tracks"]["parallel_jobs"],
        results / ".checkpoints" / "timings" / receipt_group,
        progress=progress_events(results / "logs", event_module, len(jobs), "sample"),
    )
    return (
        [path for paths, _rows in completed for path in paths],
        [row for _paths, rows in completed for row in rows],
    )


def make_c0_tracks(plan: RunPlan, results: Path, dry_run: bool = False, force: bool = False) -> None:
    """Publish all-read raw/CPM tracks immediately after final C0 BAM creation."""
    settings = plan.project["tracks"]
    logdir = results / "logs"
    enabled = bool(
        plan.project.get("modules", {}).get("tracks", True)
        and settings.get("early_c0", True)
        and settings["families"].get("all_reads", False)
    )
    if not enabled:
        event(logdir, "c0_tracks", "disabled", "Early C0 track generation is disabled")
        return
    normalizations = tuple(
        key for key in ("raw", "cpm") if settings["normalizations"].get(key, False)
    )
    if dry_run:
        event(logdir, "c0_tracks", "dry_run",
              "Would publish all-read " + "/".join(normalizations) + " strand-specific tracks")
        return
    required = ["samtools", "bedtools"]
    if settings["generate_bigwigs"]:
        required.append("bedGraphToBigWig")
    require_tools(required)
    outdir = results / "09_tracks"
    stage_receipt = outdir / ".stage_receipts" / "tracks_c0"
    signature_inputs = [
        results / "02_alignment" / sample["sample_id"] / f"{sample['sample_id']}.bam"
        for sample in plan.samples
    ]
    signature_inputs.extend(Path(plan.reference_for(sample["genome"])["chrom_sizes"]) for sample in plan.samples)
    signature = signature_for(signature_inputs, {
        "module": "tracks_c0", "samples": plan.samples,
        "normalizations": normalizations,
        "generate_bigwigs": settings["generate_bigwigs"],
        "retain_bedgraph": settings["retain_bedgraph"],
    })
    if not force and receipt_valid(stage_receipt, signature):
        event(logdir, "c0_tracks", "skipped", "Valid early C0 track receipt")
        return
    outputs, rows = _run_track_subset(
        plan, results, force, ("all_reads",), normalizations, "tracks_c0", "c0_tracks"
    )
    normalization_path = outdir / "track_normalization.c0.tsv"
    _write_normalization(normalization_path, rows)
    outputs.append(normalization_path)
    write_receipt("tracks_c0", stage_receipt, signature, outputs, ["rna-ends2tracks", "c0_tracks"])
    event(logdir, "c0_tracks", "completed",
          f"Published {len(outputs)} early C0 track and normalization deliverables")


def make_tracks(plan: RunPlan, results: Path, dry_run: bool = False, force: bool = False) -> None:
    outdir = results / "09_tracks"
    logdir = results / "logs"
    if not plan.project.get("modules", {}).get("tracks", True):
        event(logdir, "tracks", "disabled", "RUN_TRACKS=false")
        return
    if dry_run:
        enabled_families = [key for key, value in plan.project["tracks"]["families"].items() if value]
        event(logdir, "tracks", "dry_run", "Would generate: " + ", ".join(enabled_families))
        return
    settings = plan.project["tracks"]
    early_enabled = bool(settings.get("early_c0", True) and settings["families"].get("all_reads", False))
    if early_enabled:
        make_c0_tracks(plan, results, False, force)
    final_families = tuple(
        family for family, enabled in settings["families"].items()
        if enabled and not (early_enabled and family == "all_reads")
    )
    final_normalizations = tuple(
        normalization for normalization, enabled in settings["normalizations"].items() if enabled
    )
    required: list[str] = []
    if final_families:
        required.extend(["samtools", "bedtools"])
    if final_families and settings["generate_bigwigs"]:
        required.append("bedGraphToBigWig")
    final_norms = settings["normalizations"]
    needs_final_factors = (final_norms["deseq2"] or final_norms["robust_cpm"]) and (
        "filtered_ends" in final_families or "active_pas" in final_families
    )
    if needs_final_factors:
        required.append("Rscript")
    if required:
        require_tools(required)
    if needs_final_factors:
        for genome in plan.references or {plan.reference["assembly"]: plan.reference}:
            if not any(sample["genome"] == genome for sample in plan.samples):
                continue
            factor_path = results / "04_active_pas" / genome / "C4_track_size_factors.tsv"
            c4_path = results / "04_active_pas" / genome / "C4_active_pas_gene_counts.tsv"
            samples_path = results / "00_metadata" / "validated_samples.tsv"
            factor_receipt = factor_path.parent / ".track_factor_receipt"
            factor_signature = signature_for(
                [c4_path, samples_path],
                {"module": "C4_track_size_factors", "genome": genome, "estimator": "DESeq2_poscounts"},
            )
            if not receipt_valid(factor_receipt, factor_signature):
                temporary_factor = factor_path.with_name(f".{factor_path.name}.tmp")
                run([
                    "Rscript", str(workflow_asset("scripts/R/c4_track_size_factors.R")),
                    "--counts", str(c4_path), "--samples", str(samples_path),
                    "--output", str(temporary_factor),
                ], results / "logs" / "tracks" / genome / "size_factors.log", False)
                temporary_factor.replace(factor_path)
                factor_receipt.mkdir(parents=True, exist_ok=True)
                write_receipt(
                    "C4_track_size_factors", factor_receipt, factor_signature, [factor_path],
                    ["rna-ends2tracks", "C4_track_size_factors", genome],
                )
    if final_families:
        end_outputs, end_rows = _run_track_subset(
            plan, results, force, final_families, final_normalizations, "tracks_ends", "tracks"
        )
    else:
        end_outputs, end_rows = [], []
    end_normalization = outdir / "track_normalization.ends.tsv"
    _write_normalization(end_normalization, end_rows)
    end_outputs.append(end_normalization)
    end_receipt = outdir / ".stage_receipts" / "tracks_ends"
    end_signature_inputs = [
        outdir / ".receipts" / "tracks_ends" / sample["sample_id"] / "run_receipt.json"
        for sample in plan.samples
        if final_families
    ]
    end_signature = signature_for(end_signature_inputs, {
        "module": "tracks_ends", "families": final_families,
        "normalizations": final_normalizations, "samples": plan.samples,
    })
    write_receipt("tracks_ends", end_receipt, end_signature, end_outputs,
                  ["rna-ends2tracks", "tracks_ends"])
    expected = list(end_outputs)
    if needs_final_factors:
        expected.extend(
            results / "04_active_pas" / genome / "C4_track_size_factors.tsv"
            for genome in plan.references
            if any(sample["genome"] == genome for sample in plan.samples)
        )
    rows: list[dict[str, Any]] = []
    if early_enabled:
        c0_receipt = outdir / ".stage_receipts" / "tracks_c0"
        expected.extend(_receipt_outputs(c0_receipt))
        rows.extend(_read_normalization(outdir / "track_normalization.c0.tsv"))
    rows.extend(end_rows)
    normalization_path = outdir / "track_normalization.tsv"
    _write_normalization(normalization_path, rows)
    expected.append(normalization_path)
    signature_inputs = [end_receipt / "run_receipt.json"]
    if early_enabled:
        signature_inputs.append(outdir / ".stage_receipts" / "tracks_c0" / "run_receipt.json")
    signature = signature_for(signature_inputs, {"module": "tracks", "settings": plan.project["tracks"], "samples": plan.samples})
    write_receipt("tracks", outdir, signature, expected, ["rna-ends2tracks", "tracks"])
    event(logdir, "tracks", "completed", f"Published {len(expected)} track and normalization deliverables")
