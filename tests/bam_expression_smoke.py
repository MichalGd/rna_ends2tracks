"""Linux install canary: actual pysam BAM reuse and featureCounts SAF/GTF counts."""
from __future__ import annotations

import csv
import subprocess
import sys
from pathlib import Path

import pysam

from rnaends2tracks.alignment_inputs import validate_reused_alignments
from rnaends2tracks.gex_annotation import build_annotation


def check(root: Path) -> None:
    root = root.resolve()
    source = root / "source" / "fixture"
    source.mkdir(parents=True, exist_ok=True)
    bam = source / "fixture.bam"
    header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": 1000}],
              "RG": [{"ID": "fixture", "SM": "fixture"}]}
    with pysam.AlignmentFile(str(bam), "wb", header=header) as output:
        for i, position in enumerate((120, 120, 140, 280, 320, 500)):
            read = pysam.AlignedSegment()
            read.query_name = f"r{i}"
            read.flag = 0 if i == 1 else 16
            read.reference_id = 0
            read.reference_start = position - 1
            read.mapping_quality = 255
            read.cigartuples = [(0, 11), (3, 49), (0, 20)] if i == 2 else [(0, 20)]
            length = 31 if i == 2 else 20
            read.query_sequence = "A" * length
            read.query_qualities = pysam.qualitystring_to_array("I" * length)
            read.set_tag("NH", 1)
            read.set_tag("RG", "fixture")
            output.write(read)
    pysam.index(str(bam))
    sizes = root / "chrom.sizes"
    sizes.write_text("chr1\t1000\n", encoding="utf-8")
    records = validate_reused_alignments(
        {"input_mode": "bam", "alignment_source_dir": str(root / "source"), "output_dir": str(root / "results")},
        [{"sample_id": "fixture", "genome": "GRCm39"}], {"GRCm39": {"chrom_sizes": str(sizes)}},
    )
    assert records[0]["sampled_c0_records"] == 6
    gtf = root / "genes.gtf"
    features = [("A", 100, 150, "+"), ("A", 200, 300, "+"), ("B", 250, 350, "+"), ("C", 100, 200, "-")]
    gtf.write_text("".join(f'chr1\tfixture\texon\t{start}\t{end}\t.\t{strand}\t.\tgene_id "{gene}"; transcript_id "{gene}.t1"; gene_type "protein_coding";\n'
                           for gene, start, end, strand in features), encoding="utf-8")
    saf = build_annotation(gtf, sizes, root / "annotation", extension=0)["saf"]
    observed = []
    for kind, annotation, arguments in (("saf", saf, ["-F", "SAF"]), ("gtf", gtf, ["-t", "exon", "-g", "gene_id"])):
        table = root / f"counts.{kind}.tsv"
        command = ["featureCounts", "-T", "1", "-a", str(annotation), "-o", str(table),
                   "-s", "2", "--primary", "--tmpDir", str(root), *arguments, str(bam)]
        with (root / f"featureCounts.{kind}.log").open("w", encoding="utf-8") as log:
            subprocess.run(command, check=True, stdout=log, stderr=subprocess.STDOUT)
        with table.open(encoding="utf-8", newline="") as handle:
            values = {row["Geneid"]: int(row[str(bam)]) for row in csv.DictReader(
                (line for line in handle if not line.startswith("#")), delimiter="\t")}
        assert values == {"A": 2, "B": 1, "C": 1}, values
        observed.append(values)
        with Path(str(table) + ".summary").open(encoding="utf-8", newline="") as handle:
            summary = {row["Status"]: int(row[str(bam)]) for row in csv.DictReader(handle, delimiter="\t")}
        assert summary["Unassigned_Ambiguity"] == summary["Unassigned_NoFeatures"] == 1
    assert observed[0] == observed[1]
    print("PASS: actual indexed C0 BAM reuse; zero-extension SAF equals GTF; stranded counts, ambiguity and splice counting.")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise SystemExit("Usage: python tests/bam_expression_smoke.py NEW_FIXTURE_DIRECTORY")
    check(Path(sys.argv[1]))
