from __future__ import annotations

import csv
import gzip
import json
import re
import tempfile
from collections import defaultdict
from pathlib import Path
from typing import Any

from .config import RunPlan, signature_for
from .execution import run_bounded
from .external import event, progress_events, require_tools, run, run_to_path
from .preprocess import _remove_owned_temporary_tree
from .receipts import receipt_valid, write_receipt


SUMMARY_FIELDS = [
    "sample_id", "genome", "condition", "infer_undetermined_fraction",
    "infer_forward_fraction", "infer_reverse_fraction", "dominant_orientation",
    "five_prime_fraction", "three_prime_fraction", "three_to_five_ratio",
    "quantseq_three_prime_enriched", "infer_experiment_file",
    "read_distribution_file", "gene_body_coverage_file",
]


def _gtf_attributes(value: str) -> dict[str, str]:
    return {
        key: item
        for key, item in re.findall(r'([A-Za-z0-9_.-]+)\s+"([^"]*)"', value)
    }


def gtf_to_bed12(gtf: str | Path, output: Path) -> int:
    """Create a deterministic transcript BED12 reference from GTF exons."""
    source = Path(gtf)
    opener = gzip.open if source.suffix == ".gz" else open
    transcripts: dict[tuple[str, str, str], list[tuple[int, int]]] = defaultdict(list)
    names: dict[tuple[str, str, str], str] = {}
    with opener(source, "rt", encoding="utf-8") as handle:
        for number, line in enumerate(handle, 1):
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "exon":
                continue
            attributes = _gtf_attributes(fields[8])
            transcript_id = attributes.get("transcript_id", "")
            if not transcript_id:
                continue
            chrom, strand = fields[0], fields[6]
            if strand not in {"+", "-"}:
                continue
            try:
                start, end = int(fields[3]) - 1, int(fields[4])
            except ValueError as exc:
                raise RuntimeError(f"Invalid GTF coordinates at {source}:{number}") from exc
            if start < 0 or end <= start:
                raise RuntimeError(f"Invalid GTF exon interval at {source}:{number}")
            key = (chrom, strand, transcript_id)
            transcripts[key].append((start, end))
            names[key] = attributes.get("gene_name") or attributes.get("gene_id") or transcript_id

    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_name(f".{output.name}.tmp")
    written = 0
    with temporary.open("w", encoding="utf-8", newline="") as handle:
        for key in sorted(transcripts):
            chrom, strand, transcript_id = key
            exons = sorted(set(transcripts[key]))
            # BED12 blocks may not overlap. GENCODE exon records for a single
            # transcript normally do not, but fail clearly if a malformed GTF does.
            if any(left[1] > right[0] for left, right in zip(exons, exons[1:])):
                raise RuntimeError(f"Overlapping exons in transcript {transcript_id}")
            tx_start, tx_end = exons[0][0], exons[-1][1]
            label = f"{names[key]}|{transcript_id}".replace("\t", "_").replace(" ", "_")
            block_sizes = ",".join(str(end - start) for start, end in exons) + ","
            block_starts = ",".join(str(start - tx_start) for start, _end in exons) + ","
            fields = [
                chrom, str(tx_start), str(tx_end), label, "0", strand,
                str(tx_start), str(tx_end), "0", str(len(exons)), block_sizes, block_starts,
            ]
            handle.write("\t".join(fields) + "\n")
            written += 1
    if not written:
        temporary.unlink(missing_ok=True)
        raise RuntimeError(f"No transcript exons were available for BED12 generation: {source}")
    temporary.replace(output)
    return written


def _parse_infer_experiment(path: Path) -> dict[str, str]:
    result = {
        "infer_undetermined_fraction": "", "infer_forward_fraction": "",
        "infer_reverse_fraction": "", "dominant_orientation": "",
    }
    fractions: list[float] = []
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        match = re.search(r"Fraction of reads failed to determine:\s*([0-9.eE+-]+)", line)
        if match:
            result["infer_undetermined_fraction"] = match.group(1)
            continue
        match = re.search(r"Fraction of reads explained by .*:\s*([0-9.eE+-]+)", line)
        if match:
            fractions.append(float(match.group(1)))
    if fractions:
        result["infer_forward_fraction"] = f"{fractions[0]:.6f}"
    if len(fractions) > 1:
        result["infer_reverse_fraction"] = f"{fractions[1]:.6f}"
        result["dominant_orientation"] = (
            "reverse" if fractions[1] > fractions[0]
            else "forward" if fractions[0] > fractions[1]
            else "undetermined"
        )
    return result


def _gene_body_values(path: Path) -> list[float]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader, [])
        row = next((candidate for candidate in reader if candidate), [])
    if len(header) < 101 or len(row) < 101:
        raise RuntimeError(f"Invalid RSeQC gene-body coverage table: {path}")
    try:
        return [float(value) for value in row[1:101]]
    except ValueError as exc:
        raise RuntimeError(f"Non-numeric RSeQC gene-body coverage values: {path}") from exc


def _coverage_metrics(values: list[float]) -> dict[str, str]:
    total = sum(values)
    five = sum(values[:20])
    three = sum(values[-20:])
    return {
        "five_prime_fraction": f"{five / total:.6f}" if total else "",
        "three_prime_fraction": f"{three / total:.6f}" if total else "",
        "three_to_five_ratio": f"{three / five:.6f}" if five else "",
        "quantseq_three_prime_enriched": str(bool(three > five)).lower() if total else "",
    }


def _write_gene_body_outputs(
    rows: list[tuple[dict[str, str], list[float]]], table: Path, svg: Path,
) -> None:
    table.parent.mkdir(parents=True, exist_ok=True)
    with table.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["sample_id", "genome", "condition", *[f"p{i}" for i in range(1, 101)]])
        for sample, values in rows:
            maximum = max(values) if values else 0
            normalized = [value / maximum if maximum else 0 for value in values]
            writer.writerow([
                sample["sample_id"], sample.get("genome", ""), sample.get("condition", ""),
                *[f"{value:.8f}" for value in normalized],
            ])

    width, height = 1100, 620
    left, top, plot_width, plot_height = 80, 45, 760, 500
    palette = [
        "#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9",
        "#000000", "#7A5195", "#EF5675", "#2F4B7C", "#59A14F", "#F28E2B",
    ]
    elements = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        f'<line x1="{left}" y1="{top + plot_height}" x2="{left + plot_width}" y2="{top + plot_height}" stroke="#333"/>',
        f'<line x1="{left}" y1="{top}" x2="{left}" y2="{top + plot_height}" stroke="#333"/>',
        f'<text x="{left + plot_width / 2}" y="{height - 22}" text-anchor="middle" font-family="sans-serif">Gene-body percentile (5-prime to 3-prime)</text>',
        f'<text x="20" y="{top + plot_height / 2}" transform="rotate(-90 20 {top + plot_height / 2})" text-anchor="middle" font-family="sans-serif">Relative coverage</text>',
        f'<text x="{left}" y="{top + plot_height + 22}" text-anchor="middle" font-family="sans-serif">1</text>',
        f'<text x="{left + plot_width}" y="{top + plot_height + 22}" text-anchor="middle" font-family="sans-serif">100</text>',
    ]
    for index, (sample, values) in enumerate(rows):
        maximum = max(values) if values else 0
        normalized = [value / maximum if maximum else 0 for value in values]
        points = " ".join(
            f"{left + plot_width * position / 99:.2f},{top + plot_height * (1 - value):.2f}"
            for position, value in enumerate(normalized)
        )
        color = palette[index % len(palette)]
        elements.append(f'<polyline points="{points}" fill="none" stroke="{color}" stroke-width="1.6"/>')
        legend_y = top + 18 + index * 25
        elements.extend([
            f'<line x1="870" y1="{legend_y - 5}" x2="895" y2="{legend_y - 5}" stroke="{color}" stroke-width="3"/>',
            f'<text x="902" y="{legend_y}" font-size="12" font-family="sans-serif">{sample["sample_id"]}</text>',
        ])
    elements.append("</svg>")
    svg.write_text("\n".join(elements) + "\n", encoding="utf-8")


def rseqc(plan: RunPlan, results: Path, dry_run: bool = False, force: bool = False) -> None:
    settings = plan.project.get("rseqc", {})
    module_dir = results / "01_qc" / "rseqc"
    log_dir = results / "logs"
    if not settings.get("enabled", True):
        event(log_dir, "rseqc", "disabled", "RUN_RSEQC=false")
        return
    if dry_run:
        event(
            log_dir, "rseqc", "dry_run",
            f"Would run annotation-aware RSeQC for {len(plan.samples)} samples",
        )
        return

    enabled = [
        name for name, value in (
            ("infer_experiment", settings.get("infer_experiment", True)),
            ("read_distribution", settings.get("read_distribution", True)),
            ("gene_body_coverage", settings.get("gene_body_coverage", True)),
        ) if value
    ]
    if not enabled:
        raise RuntimeError("RUN_RSEQC=true requires at least one enabled RSeQC analysis")
    tools = []
    if "infer_experiment" in enabled:
        tools.append("infer_experiment.py")
    if "read_distribution" in enabled:
        tools.append("read_distribution.py")
    if "gene_body_coverage" in enabled:
        tools.append("geneBody_coverage.py")
    if settings.get("multiqc", True):
        tools.append("multiqc")
    require_tools(tools)

    bams = [results / "02_alignment" / sample["sample_id"] / f"{sample['sample_id']}.bam" for sample in plan.samples]
    references: dict[str, Path] = {}
    signature_inputs: list[Path] = list(bams)
    for genome in sorted(plan.references):
        reference = plan.reference_for(genome)
        configured = str(reference.get("rseqc_bed", "") or "")
        if configured:
            bed = Path(configured)
            signature_inputs.append(bed)
        else:
            bed = module_dir / "references" / f"{genome}.{reference.get('release', 'annotation')}.bed12"
            signature_inputs.append(Path(reference["gtf"]))
        references[genome] = bed

    signature = signature_for(signature_inputs, {
        "module": "rseqc", "settings": settings,
        "samples": [{key: sample.get(key, "") for key in ("sample_id", "genome", "condition")} for sample in plan.samples],
    })
    summary_path = module_dir / "rseqc_summary.tsv"
    coverage_table = module_dir / "gene_body_coverage.tsv"
    coverage_svg = module_dir / "gene_body_coverage.svg"
    multiqc_report = module_dir / "multiqc" / "multiqc_report.html"
    expected = [summary_path]
    expected.extend(
        bed for genome, bed in references.items()
        if not str(plan.reference_for(genome).get("rseqc_bed", "") or "")
    )
    for sample in plan.samples:
        sample_id = sample["sample_id"]
        if "infer_experiment" in enabled:
            expected.append(module_dir / "infer_experiment" / f"{sample_id}.infer_experiment.txt")
        if "read_distribution" in enabled:
            expected.append(module_dir / "read_distribution" / f"{sample_id}.read_distribution.txt")
        if "gene_body_coverage" in enabled:
            expected.append(module_dir / "gene_body" / f"{sample_id}.geneBodyCoverage.txt")
    if "gene_body_coverage" in enabled:
        expected.extend([coverage_table, coverage_svg])
    if settings.get("multiqc", True):
        expected.append(multiqc_report)
    if not force and receipt_valid(module_dir, signature):
        event(log_dir, "rseqc", "skipped", "Valid matching RSeQC receipt")
        return

    for genome, bed in references.items():
        if not bed.is_file():
            count = gtf_to_bed12(plan.reference_for(genome)["gtf"], bed)
            event(log_dir, "rseqc", "progress", f"Generated {genome} BED12 reference with {count} transcripts")

    jobs: list[tuple[str, Any]] = []
    for sample in plan.samples:
        sample_id = sample["sample_id"]
        bam = results / "02_alignment" / sample_id / f"{sample_id}.bam"
        bed = references[sample["genome"]]

        def worker(sample=sample, sample_id=sample_id, bam=bam, bed=bed):
            sample_log = log_dir / "rseqc" / f"{sample_id}.log"
            row = {field: "" for field in SUMMARY_FIELDS}
            row.update({
                "sample_id": sample_id, "genome": sample.get("genome", ""),
                "condition": sample.get("condition", ""),
            })
            if "infer_experiment" in enabled:
                target = module_dir / "infer_experiment" / f"{sample_id}.infer_experiment.txt"
                run_to_path([
                    "infer_experiment.py", "-i", str(bam), "-r", str(bed),
                    "-s", str(settings.get("sample_reads", 200000)),
                ], target, sample_log)
                row.update(_parse_infer_experiment(target))
                row["infer_experiment_file"] = str(target)
            if "read_distribution" in enabled:
                target = module_dir / "read_distribution" / f"{sample_id}.read_distribution.txt"
                run_to_path(["read_distribution.py", "-i", str(bam), "-r", str(bed)], target, sample_log)
                row["read_distribution_file"] = str(target)
            values: list[float] = []
            if "gene_body_coverage" in enabled:
                prefix = module_dir / "gene_body" / sample_id
                prefix.parent.mkdir(parents=True, exist_ok=True)
                run([
                    "geneBody_coverage.py", "-i", str(bam), "-r", str(bed),
                    "-l", str(settings.get("minimum_transcript_length", 100)),
                    "-f", "png", "-o", str(prefix),
                ], sample_log)
                target = prefix.with_name(prefix.name + ".geneBodyCoverage.txt")
                if not target.is_file():
                    raise RuntimeError(f"RSeQC did not create gene-body coverage output for {sample_id}: {target}")
                values = _gene_body_values(target)
                row.update(_coverage_metrics(values))
                row["gene_body_coverage_file"] = str(target)
            return row, values

        jobs.append((sample_id, worker))

    results_by_sample = run_bounded(
        "rseqc_samples", jobs, plan.project["resources"]["rseqc"]["parallel_jobs"],
        results / ".checkpoints" / "timings" / "rseqc",
        progress=progress_events(log_dir, "rseqc", len(jobs), "sample"),
    )
    rows = [row for row, _values in results_by_sample]
    with summary_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=SUMMARY_FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    if "gene_body_coverage" in enabled:
        _write_gene_body_outputs(
            [(sample, values) for sample, (_row, values) in zip(plan.samples, results_by_sample)],
            coverage_table, coverage_svg,
        )

    if settings.get("multiqc", True):
        temp_parent = results / ".checkpoints" / "multiqc_tmp"
        temp_parent.mkdir(parents=True, exist_ok=True)
        temporary = Path(tempfile.mkdtemp(prefix="rseqc.", dir=temp_parent))
        try:
            multiqc_sources = [
                module_dir / directory for analysis, directory in (
                    ("infer_experiment", "infer_experiment"),
                    ("read_distribution", "read_distribution"),
                    ("gene_body_coverage", "gene_body"),
                ) if analysis in enabled
            ]
            run([
                "multiqc", "--force", "--no-clean-up", "--outdir", str(module_dir / "multiqc"),
                *map(str, multiqc_sources),
            ], log_dir / "rseqc" / "multiqc.log", env={"TMPDIR": str(temporary)})
        finally:
            _remove_owned_temporary_tree(temporary)

    write_receipt("rseqc", module_dir, signature, expected, ["rna-ends2tracks", "rseqc"])
    event(
        log_dir, "rseqc", "completed",
        f"RSeQC completed for {len(rows)} samples; analyses={','.join(enabled)}",
    )
