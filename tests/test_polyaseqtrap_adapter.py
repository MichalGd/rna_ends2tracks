from __future__ import annotations

import csv
import hashlib
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import pysam

from rnaends2tracks.polyaseqtrap_adapter import (
    Candidate,
    cluster_candidates,
    extract_end_counts,
    parser,
    parse_deepip,
    verify_installation,
    write_outputs,
)


class PolyAseqTrapAdapterTests(unittest.TestCase):
    def test_pilot_mode_does_not_require_an_already_accepted_manifest(self) -> None:
        args = parser().parse_args([
            "--bam-manifest", "bams.tsv", "--fasta", "genome.fa", "--gtf", "genes.gtf",
            "--species", "mouse", "--assembly", "GRCm39", "--outdir", "pilot", "--pilot-mode",
        ])
        self.assertTrue(args.pilot_mode)
        self.assertEqual(args.validation_manifest, "")

    def test_pilot_installation_verification_checks_pinned_assets(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            model = root / "mouse.hdf5"; model.write_bytes(b"model")
            script = root / "DeepIP_test.py"; script.write_text("# pilot\n", encoding="utf-8")
            digest = hashlib.sha256(model.read_bytes()).hexdigest()
            installation = {
                "engine": {"source_commit": "176ea2884ff1c6be7c64bc44fa7661d82d90e718"},
                "deepip": {"source_commit": "988564875d002b6d5d48d8dfb228cba3492dd776",
                           "script": str(script)},
                "models": {"mouse": {"path": str(model), "sha256": digest}},
                "environment": {"sha256": "a" * 64},
            }
            with patch("rnaends2tracks.polyaseqtrap_adapter.DEEPIP_MODEL_SHA256",
                       {"human": "b" * 64, "mouse": digest}):
                observed_model, observed_script, environment = verify_installation(installation, "mouse")
            self.assertEqual(observed_model, model)
            self.assertEqual(observed_script, script)
            self.assertEqual(environment, "a" * 64)

    def test_project_clustering_is_weighted_and_strand_deterministic(self) -> None:
        rows = [
            Candidate("chr1", "+", 100, {"A": 5}, {"L1"}, False, 1),
            Candidate("chr1", "+", 110, {"B": 5}, {"L2"}, True, 1),
            Candidate("chr1", "-", 200, {"A": 4}, {"L1"}, False, 1),
            Candidate("chr1", "-", 210, {"B": 4}, {"L2"}, False, 1),
        ]
        result = cluster_candidates(rows, 24)
        self.assertEqual([row.position for row in result], [110, 200])
        self.assertEqual(result[0].sample_counts, {"A": 5, "B": 5})
        self.assertTrue(result[0].repeat_detected)

    def test_endpoint_extraction_retains_duplicate_flagged_records(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            bam = root / "input.bam"
            header = {"HD": {"VN": "1.6"}, "SQ": [{"SN": "chr1", "LN": 1000}]}
            with pysam.AlignmentFile(bam, "wb", header=header) as handle:
                for flag, start in ((0, 100), (1024, 100), (16, 200)):
                    read = pysam.AlignedSegment()
                    read.query_name = f"read_{flag}_{start}"
                    read.query_sequence = "A" * 50
                    read.flag = flag
                    read.reference_id = 0
                    read.reference_start = start
                    read.mapping_quality = 255
                    read.cigar = ((0, 50),)
                    read.query_qualities = pysam.qualitystring_to_array("I" * 50)
                    handle.write(read)
            output = root / "ends.tsv"
            audit = extract_end_counts(bam, output)
            self.assertEqual(audit["records_written"], 3)
            self.assertEqual(audit["duplicate_flagged_records_retained"], 1)
            with output.open(encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(sum(int(row["count"]) for row in rows), 3)
            self.assertIn({"chrom": "chr1", "position": "100", "strand": "-", "count": "2"}, rows)
            self.assertIn({"chrom": "chr1", "position": "249", "strand": "+", "count": "1"}, rows)

    def test_deepip_csv_parsing_uses_predict_label_and_strips_fasta_label(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "predictions.csv"
            path.write_text(
                "title,score,true_label,predict_label,res\n"
                ">candidate_000000001/0,0.1,0,0,TN\n"
                ">candidate_000000002/1,0.9,1,1,TP\n",
                encoding="utf-8",
            )
            result = parse_deepip(path, {"candidate_000000001", "candidate_000000002"})
            self.assertEqual(result["candidate_000000001"][0], 0)
            self.assertEqual(result["candidate_000000002"][0], 1)

    def test_intragenic_pas_is_retained_and_annotated(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            gtf = root / "genes.gtf"
            gtf.write_text(
                'chr1\ttest\tgene\t101\t900\t.\t+\t.\tgene_id "gene1";\n'
                'chr1\ttest\texon\t101\t200\t.\t+\t.\tgene_id "gene1";\n'
                'chr1\ttest\texon\t801\t900\t.\t+\t.\tgene_id "gene1";\n',
                encoding="utf-8",
            )
            candidate = Candidate("chr1", "+", 500, {"A": 8, "B": 7}, {"QuantSeq"}, False, 2)
            write_outputs([candidate], ["A", "B"], gtf, root / "out", [], {"assembly": "GRCh38"})
            with (root / "out" / "pas_catalog.tsv").open(encoding="utf-8") as handle:
                row = next(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(row["gene_id"], "gene1")
            self.assertEqual(row["feature_class"], "intron")
            self.assertEqual(row["start"], "500")


if __name__ == "__main__":
    unittest.main()
