import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from rnaends2tracks.config import RunPlan
from rnaends2tracks.rseqc import (
    _coverage_metrics, _gene_body_values, _parse_infer_experiment,
    gtf_to_bed12, rseqc,
)


class RSeQCTests(unittest.TestCase):
    def test_gtf_to_bed12_preserves_transcript_blocks_and_strand(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            gtf = root / "annotation.gtf"
            gtf.write_text(
                'chr1\tsrc\texon\t101\t150\t.\t+\t.\tgene_id "g1"; transcript_id "t1"; gene_name "Gene1";\n'
                'chr1\tsrc\texon\t201\t250\t.\t+\t.\tgene_id "g1"; transcript_id "t1"; gene_name "Gene1";\n'
                'chr2\tsrc\texon\t401\t430\t.\t-\t.\tgene_id "g2"; transcript_id "t2";\n',
                encoding="utf-8",
            )
            bed = root / "annotation.bed12"
            self.assertEqual(gtf_to_bed12(gtf, bed), 2)
            rows = [line.split("\t") for line in bed.read_text(encoding="utf-8").splitlines()]
            self.assertEqual(rows[0][:6], ["chr1", "100", "250", "Gene1|t1", "0", "+"])
            self.assertEqual(rows[0][9:], ["2", "50,50,", "0,100,"])
            self.assertEqual(rows[1][5], "-")

    def test_rseqc_parsers_report_orientation_and_three_prime_enrichment(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            inferred = root / "infer.txt"
            inferred.write_text(
                "Fraction of reads failed to determine: 0.02\n"
                'Fraction of reads explained by "1++,1--,2+-,2-+": 0.10\n'
                'Fraction of reads explained by "1+-,1-+,2++,2--": 0.88\n',
                encoding="utf-8",
            )
            parsed = _parse_infer_experiment(inferred)
            self.assertEqual(parsed["dominant_orientation"], "reverse")
            self.assertEqual(parsed["infer_reverse_fraction"], "0.880000")
            coverage = root / "coverage.txt"
            coverage.write_text(
                "sample\t" + "\t".join(str(index) for index in range(1, 101)) + "\n"
                "S1\t" + "\t".join(str(index) for index in range(1, 101)) + "\n",
                encoding="utf-8",
            )
            values = _gene_body_values(coverage)
            metrics = _coverage_metrics(values)
            self.assertEqual(len(values), 100)
            self.assertEqual(metrics["quantseq_three_prime_enriched"], "true")
            self.assertGreater(float(metrics["three_to_five_ratio"]), 1)

    def test_stage_generates_summary_plot_multiqc_and_receipt(self):
        with tempfile.TemporaryDirectory() as temporary:
            results = Path(temporary) / "results"
            gtf = Path(temporary) / "annotation.gtf"
            gtf.write_text(
                'chr1\tsrc\texon\t101\t500\t.\t+\t.\tgene_id "g1"; transcript_id "t1";\n',
                encoding="utf-8",
            )
            samples = [
                {"sample_id": "S1", "genome": "GRCm39", "condition": "control"},
                {"sample_id": "S2", "genome": "GRCm39", "condition": "treated"},
            ]
            for sample in samples:
                bam = results / "02_alignment" / sample["sample_id"] / f"{sample['sample_id']}.bam"
                bam.parent.mkdir(parents=True, exist_ok=True)
                bam.write_bytes(b"bam")
            reference = {"assembly": "GRCm39", "release": "vM31", "gtf": str(gtf), "rseqc_bed": ""}
            project = {
                "rseqc": {
                    "enabled": True, "infer_experiment": True, "read_distribution": True,
                    "gene_body_coverage": True, "multiqc": True,
                    "sample_reads": 1000, "minimum_transcript_length": 100,
                },
                "resources": {"rseqc": {"parallel_jobs": 2}},
            }
            plan = RunPlan(project, samples, samples, [], reference, {"GRCm39": reference})

            def fake_run_to_path(command, target, log, **_kwargs):
                target.parent.mkdir(parents=True, exist_ok=True)
                if command[0] == "infer_experiment.py":
                    target.write_text(
                        "Fraction of reads failed to determine: 0.01\n"
                        'Fraction of reads explained by "forward": 0.05\n'
                        'Fraction of reads explained by "reverse": 0.94\n',
                        encoding="utf-8",
                    )
                else:
                    target.write_text("Total Reads\t1000\n", encoding="utf-8")

            def fake_run(command, log, **_kwargs):
                if command[0] == "geneBody_coverage.py":
                    prefix = Path(command[command.index("-o") + 1])
                    target = prefix.with_name(prefix.name + ".geneBodyCoverage.txt")
                    target.write_text(
                        "sample\t" + "\t".join(str(index) for index in range(1, 101)) + "\n"
                        + prefix.name + "\t" + "\t".join(str(index) for index in range(1, 101)) + "\n",
                        encoding="utf-8",
                    )
                elif command[0] == "multiqc":
                    outdir = Path(command[command.index("--outdir") + 1])
                    outdir.mkdir(parents=True, exist_ok=True)
                    (outdir / "multiqc_report.html").write_text("<html></html>\n", encoding="utf-8")

            with (
                patch("rnaends2tracks.rseqc.require_tools"),
                patch("rnaends2tracks.rseqc.run_to_path", side_effect=fake_run_to_path),
                patch("rnaends2tracks.rseqc.run", side_effect=fake_run),
            ):
                rseqc(plan, results)

            root = results / "01_qc" / "rseqc"
            self.assertTrue((root / "rseqc_summary.tsv").is_file())
            self.assertTrue((root / "gene_body_coverage.svg").is_file())
            self.assertTrue((root / "multiqc" / "multiqc_report.html").is_file())
            self.assertTrue((root / "run_receipt.json").is_file())
            self.assertIn("RSeQC completed for 2 samples", (results / "rna_ends2tracks.log").read_text())


if __name__ == "__main__":
    unittest.main()
