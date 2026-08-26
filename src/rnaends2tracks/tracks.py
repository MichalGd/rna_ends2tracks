from __future__ import annotations

import csv
import gzip
import json
import subprocess
from pathlib import Path
from typing import Any

from .config import RunPlan, signature_for
from .execution import run_bounded
from .external import event, require_tools, run
from .paths import workflow_asset
from .receipts import receipt_valid, write_receipt

END_FAMILIES = {
    "exact_ends": "C1_exact_ends.tsv.gz",
    "filtered_ends": "C2_filtered_ends.tsv.gz",
    "rejected_ends": "C2R_internal_priming_rejects.tsv.gz",
}


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
    completed = subprocess.run(
        ["samtools", "view", "-@", str(threads), "-c", str(bam)], check=True,
        capture_output=True, text=True,
    )
    log.parent.mkdir(parents=True, exist_ok=True)
    log.write_text(completed.stderr, encoding="utf-8")
    return int(completed.stdout.strip())


def _all_read_bedgraph(bam: Path, strand: str, output: Path, scale: float, threads: int, log: Path) -> Path:
    strand_bam = output.with_suffix(".strand.bam")
    flag_args = ["-f", "16"] if strand == "plus" else ["-F", "16"]
    run(["samtools", "view", "-@", str(threads), "-b", *flag_args, "-o", str(strand_bam), str(bam)], log, False)
    positive = output if strand == "plus" else output.with_suffix(".positive.bedGraph")
    with positive.open("w", encoding="utf-8") as handle:
        subprocess.run(
            ["bedtools", "genomecov", "-ibam", str(strand_bam), "-bg", "-scale", f"{scale:.15g}"],
            stdout=handle, stderr=subprocess.PIPE, check=True, text=True,
        )
    if strand == "minus":
        _negated(positive, output)
    return strand_bam


def _sample_tracks(plan: RunPlan, results: Path, sample: dict[str, str], force: bool) -> tuple[list[Path], list[dict[str, Any]]]:
    sample_id, genome = sample["sample_id"], sample["genome"]
    reference = plan.reference_for(genome)
    settings = plan.project["tracks"]
    resource = plan.project["resources"]["tracks"]
    track_root = results / "09_tracks"
    intermediate = track_root / ".intermediate" / sample_id
    receipt_dir = track_root / ".receipts" / sample_id
    bam = results / "02_alignment" / sample_id / f"{sample_id}.bam"
    size_factor_path = results / "04_active_pas" / genome / "C4_track_size_factors.tsv"
    inputs: list[Path] = [bam, Path(reference["chrom_sizes"])]
    needs_factors = (
        settings["normalizations"]["deseq2"] or settings["normalizations"]["robust_cpm"]
    ) and (settings["families"]["active_pas"] or settings["families"]["filtered_ends"])
    if needs_factors:
        inputs.append(size_factor_path)
    for family, filename in END_FAMILIES.items():
        if settings["families"][family]:
            inputs.append(results / "03_exact_ends" / genome / sample_id / filename)
    if settings["families"]["active_pas"]:
        inputs.extend([
            results / "04_active_pas" / genome / "active_pas_catalog.tsv",
            results / "04_active_pas" / genome / "C3_active_pas_counts.tsv",
        ])
    signature = signature_for(inputs, {"module": "tracks", "sample": sample, "settings": settings})
    if not force and receipt_valid(receipt_dir, signature):
        receipt = receipt_dir / "run_receipt.json"
        fragment = receipt_dir / "normalization.tsv"
        receipt_payload = json.loads(receipt.read_text(encoding="utf-8"))
        with fragment.open(encoding="utf-8", newline="") as handle:
            return [Path(row["path"]) for row in receipt_payload["outputs"]], list(csv.DictReader(handle, delimiter="\t"))
    factors = _size_factors(size_factor_path)[sample_id] if size_factor_path.is_file() else {}
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
    if settings["families"]["all_reads"]:
        denominators["all_reads"] = _bam_count(bam, results / "logs" / "tracks" / f"{sample_id}.count.log", resource["samtools_threads"])
    for family, filename in END_FAMILIES.items():
        if settings["families"][family]:
            rows = _end_counts(results / "03_exact_ends" / genome / sample_id / filename)
            family_rows[family] = rows
            denominators[family] = sum(int(row["count"]) for row in rows)
    if settings["families"]["rejected_ends"]:
        if "exact_ends" not in denominators:
            c1_rows = _end_counts(results / "03_exact_ends" / genome / sample_id / END_FAMILIES["exact_ends"])
            denominators["exact_ends"] = sum(int(row["count"]) for row in c1_rows)
        denominators["rejected_ends"] = denominators["exact_ends"]
    if settings["families"]["active_pas"]:
        rows, assigned = _active_counts(
            results / "04_active_pas" / genome / "active_pas_catalog.tsv",
            results / "04_active_pas" / genome / "C3_active_pas_counts.tsv", sample_id,
        )
        family_rows["active_pas"] = rows
        denominators["active_pas"] = assigned

    for family, enabled in settings["families"].items():
        if not enabled:
            continue
        normalizations: list[tuple[str, float]] = []
        if settings["normalizations"]["raw"]:
            normalizations.append(("raw", 1.0))
        if settings["normalizations"]["cpm"]:
            denominator = denominators.get(family, 0)
            normalizations.append(("cpm", 1_000_000.0 / denominator if denominator else 0.0))
        if family in {"filtered_ends", "active_pas"}:
            if settings["normalizations"]["deseq2"]:
                normalizations.append(("deseq2", factors["deseq2"]))
            if settings["normalizations"]["robust_cpm"]:
                normalizations.append(("robust_cpm", factors["robust_cpm"]))
        for normalization, scale in normalizations:
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
                        bam, strand, bedgraph, scale, resource["samtools_threads"],
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
        fields = ["sample_id", "genome", "family", "normalization", "denominator", "scale", "count_universe"]
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(normalization_rows)
    expected.append(fragment)
    write_receipt("tracks_sample", receipt_dir, signature, expected, ["rna-ends2tracks", "tracks", sample_id])
    return expected, normalization_rows


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
    required = ["samtools", "bedtools"]
    if plan.project["tracks"]["generate_bigwigs"]:
        required.append("bedGraphToBigWig")
    final_norms = plan.project["tracks"]["normalizations"]
    needs_final_factors = (final_norms["deseq2"] or final_norms["robust_cpm"]) and (
        plan.project["tracks"]["families"]["filtered_ends"] or
        plan.project["tracks"]["families"]["active_pas"]
    )
    if needs_final_factors:
        required.append("Rscript")
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
    jobs = [(sample["sample_id"], lambda sample=sample: _sample_tracks(plan, results, sample, force)) for sample in plan.samples]
    completed = run_bounded(
        "tracks", jobs, plan.project["resources"]["tracks"]["parallel_jobs"],
        results / ".checkpoints" / "timings" / "tracks",
    )
    expected = [path for paths, _ in completed for path in paths]
    if needs_final_factors:
        expected.extend(
            results / "04_active_pas" / genome / "C4_track_size_factors.tsv"
            for genome in plan.references
            if any(sample["genome"] == genome for sample in plan.samples)
        )
    rows = [row for _, normalization in completed for row in normalization]
    normalization_path = outdir / "track_normalization.tsv"
    with normalization_path.open("w", encoding="utf-8", newline="") as handle:
        fields = ["sample_id", "genome", "family", "normalization", "denominator", "scale", "count_universe"]
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)
    expected.append(normalization_path)
    signature_inputs = [results / "02_alignment" / sample["sample_id"] / f"{sample['sample_id']}.bam" for sample in plan.samples]
    signature = signature_for(signature_inputs, {"module": "tracks", "settings": plan.project["tracks"], "samples": plan.samples})
    write_receipt("tracks", outdir, signature, expected, ["rna-ends2tracks", "tracks"])
    event(logdir, "tracks", "completed", f"Published {len(expected)} track and normalization deliverables")
