from __future__ import annotations

import subprocess
from pathlib import Path

from .config import RunPlan, signature_for
from .external import event, require_tools, run
from .receipts import receipt_valid, write_receipt


def _negate_bedgraph(source: Path, target: Path) -> None:
    with source.open(encoding="utf-8") as inp, target.open("w", encoding="utf-8") as out:
        for line in inp:
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 4:
                fields[3] = str(-float(fields[3]))
            out.write("\t".join(fields) + "\n")


def make_tracks(plan: RunPlan, results: Path, dry_run: bool = False, force: bool = False) -> None:
    outdir = results / "07_tracks"
    logdir = results / "provenance" / "logs" / "tracks"
    outdir.mkdir(parents=True, exist_ok=True)
    bams = [results / "02_alignment" / sample["sample_id"] / f"{sample['sample_id']}.bam" for sample in plan.samples]
    expected = [outdir / sample["sample_id"] / f"{sample['sample_id']}.{strand}.bw"
                for sample in plan.samples for strand in ("transcript_plus", "transcript_minus")]
    signature = signature_for([*bams, plan.reference["chrom_sizes"]], {"module": "tracks", "samples": plan.samples}) if not dry_run else "dry-run"
    if not force and not dry_run and receipt_valid(outdir, signature):
        event(results / "provenance" / "logs", "tracks", "skipped", "Valid matching receipt")
        return
    if not dry_run:
        require_tools(["samtools", "bedtools", "bedGraphToBigWig"])
    threads = str(max(1, int(plan.project.get("resources", {}).get("threads", 8))))
    for sample in plan.samples:
        sample_id = sample["sample_id"]
        bam = results / "02_alignment" / sample_id / f"{sample_id}.bam"
        sample_dir = outdir / sample_id
        sample_dir.mkdir(parents=True, exist_ok=True)
        plus_bam = sample_dir / "transcript_plus.bam"
        minus_bam = sample_dir / "transcript_minus.bam"
        plus_bg = sample_dir / "transcript_plus.bedGraph"
        minus_raw = sample_dir / "transcript_minus.positive.bedGraph"
        minus_bg = sample_dir / "transcript_minus.bedGraph"
        log = logdir / f"{sample_id}.log"
        # REV R1 aligns opposite the transcript: reverse BAM alignments are plus-transcript reads.
        run(["samtools", "view", "-@", threads, "-b", "-f", "16", "-F", "2308", "-o", str(plus_bam), str(bam)], log, dry_run)
        run(["samtools", "view", "-@", threads, "-b", "-F", "2324", "-o", str(minus_bam), str(bam)], log, dry_run)
        if dry_run:
            run(["bedtools", "genomecov", "-ibam", str(plus_bam), "-bg"], log, True)
            run(["bedtools", "genomecov", "-ibam", str(minus_bam), "-bg"], log, True)
        else:
            with plus_bg.open("w", encoding="utf-8") as handle:
                subprocess.run(["bedtools", "genomecov", "-ibam", str(plus_bam), "-bg"], stdout=handle, check=True, text=True)
            with minus_raw.open("w", encoding="utf-8") as handle:
                subprocess.run(["bedtools", "genomecov", "-ibam", str(minus_bam), "-bg"], stdout=handle, check=True, text=True)
            _negate_bedgraph(minus_raw, minus_bg)
        run(["bedGraphToBigWig", str(plus_bg), plan.reference["chrom_sizes"], str(sample_dir / f"{sample_id}.transcript_plus.bw")], log, dry_run)
        run(["bedGraphToBigWig", str(minus_bg), plan.reference["chrom_sizes"], str(sample_dir / f"{sample_id}.transcript_minus.bw")], log, dry_run)
    if not dry_run:
        write_receipt("tracks", outdir, signature, expected, ["rna-ends2tracks", "tracks"])
    event(results / "provenance" / "logs", "tracks", "dry_run" if dry_run else "completed", "REV transcript-strand browser tracks")
