from __future__ import annotations

import csv
import gzip
import hashlib
import json
import os
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import pysam

from rnaends2tracks.polyaseqtrap_adapter import (
    Candidate,
    cluster_candidates,
    extract_end_counts,
    environment_executable,
    parser,
    parse_deepip,
    r_subprocess_environment,
    reuse_exact_end_counts,
    verify_installation,
    write_outputs,
)


class PolyAseqTrapAdapterTests(unittest.TestCase):
    def test_receipt_validated_c1_plus_c1s_reconstructs_raw_endpoint_universe(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            c1 = root / "C1_exact_ends.tsv.gz"
            c1s = root / "C1S_uncertain_ends.tsv.gz"
            audit = root / "end_audit.json"
            for path, rows in (
                (c1, ["chr1\t100\t101\t+\t3", "chr1\t200\t201\t-\t2"]),
                (c1s, ["chr1\t100\t101\t+\t1", "chr2\t300\t301\t+\t4"]),
            ):
                with gzip.open(path, "wt", encoding="utf-8") as handle:
                    handle.write("chrom\tstart\tend\tstrand\tcount\n" + "\n".join(rows) + "\n")
            audit.write_text(json.dumps({
                "sample_id": "S1", "C0": 10, "C1": 5, "C1S": 5,
                "duplicate_flagged": 4, "duplicate_flagged_C0": 2,
            }), encoding="utf-8")
            receipt_dir = root / ".receipt"
            receipt_dir.mkdir()
            outputs = []
            for path in (c1, c1s, audit):
                outputs.append({
                    "path": str(path.resolve()), "size": path.stat().st_size,
                    "validation": "sha256", "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
                })
            receipt = receipt_dir / "run_receipt.json"
            receipt.write_text(json.dumps({
                "module": "exact_ends_sample", "exit_status": 0, "outputs": outputs,
            }), encoding="utf-8")
            target = root / "adapter_ends.tsv"
            observed = reuse_exact_end_counts({
                "sample_id": "S1", "c1": str(c1), "c1s": str(c1s),
                "exact_end_audit": str(audit), "exact_end_receipt": str(receipt),
            }, target)
            self.assertEqual(observed["records_written"], 10)
            self.assertEqual(observed["end_soft_clipped_records_included"], 5)
            self.assertEqual(observed["duplicate_flagged_records_retained"], 2)
            self.assertEqual(observed["duplicate_flagged_count_scope"], "C0")
            self.assertEqual(observed["endpoint_source"], "receipt_validated_C1_plus_C1S")
            with target.open(encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertIn({"chrom": "chr1", "position": "100", "strand": "+", "count": "4"}, rows)
            self.assertEqual(sum(int(row["count"]) for row in rows), 10)

    def test_reusable_endpoint_checksum_mismatch_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            paths = [root / name for name in ("C1_exact_ends.tsv.gz", "C1S_uncertain_ends.tsv.gz")]
            for path in paths:
                with gzip.open(path, "wt", encoding="utf-8") as handle:
                    handle.write("chrom\tstart\tend\tstrand\tcount\nchr1\t1\t2\t+\t1\n")
            audit = root / "end_audit.json"
            audit.write_text('{"sample_id":"S1","C0":2,"C1":1,"C1S":1}', encoding="utf-8")
            receipt_dir = root / ".receipt"; receipt_dir.mkdir()
            outputs = [{
                "path": str(path.resolve()), "size": path.stat().st_size,
                "validation": "sha256", "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            } for path in (*paths, audit)]
            receipt = receipt_dir / "run_receipt.json"
            receipt.write_text(json.dumps({
                "module": "exact_ends_sample", "exit_status": 0, "outputs": outputs,
            }), encoding="utf-8")
            audit.write_text('{"sample_id":"S1","C0":3,"C1":1,"C1S":2}', encoding="utf-8")
            with self.assertRaisesRegex(RuntimeError, "checksum"):
                reuse_exact_end_counts({
                    "sample_id": "S1", "c1": str(paths[0]), "c1s": str(paths[1]),
                    "exact_end_audit": str(audit), "exact_end_receipt": str(receipt),
                }, root / "out.tsv")

    def test_rscript_is_resolved_beside_running_apa_b_python(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            python = root / "python"
            rscript = root / "Rscript"
            python.write_text("", encoding="utf-8")
            rscript.write_text("", encoding="utf-8")
            with patch("rnaends2tracks.polyaseqtrap_adapter.sys.executable", str(python)):
                self.assertEqual(environment_executable("Rscript"), str(rscript.resolve()))

    def test_r_environment_removes_parent_overrides_and_prioritizes_apa_b_bin(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            python = Path(directory) / "python"
            python.write_text("", encoding="utf-8")
            parent = {
                "PATH": os.pathsep.join(("parent", "bin")),
                "R_HOME": "wrong", "R_LIBS": "wrong", "R_LIBS_USER": "wrong",
            }
            with patch("rnaends2tracks.polyaseqtrap_adapter.sys.executable", str(python)), \
                    patch.dict("rnaends2tracks.polyaseqtrap_adapter.os.environ", parent, clear=True):
                environment = r_subprocess_environment()
            self.assertEqual(environment["PATH"].split(os.pathsep)[0], str(python.parent.resolve()))
            self.assertNotIn("R_HOME", environment)
            self.assertNotIn("R_LIBS", environment)
            self.assertNotIn("R_LIBS_USER", environment)
            self.assertEqual(environment["OMP_NUM_THREADS"], "1")
            self.assertEqual(environment["OPENBLAS_NUM_THREADS"], "1")

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
                "workflow_adapter": {"release": "v0.1.0-alpha.10.post4", "source_commit": "7654321"},
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

    def test_pilot_installation_rejects_missing_workflow_adapter_pin(self) -> None:
        with self.assertRaisesRegex(RuntimeError, "workflow-adapter"):
            verify_installation({
                "engine": {"source_commit": "176ea2884ff1c6be7c64bc44fa7661d82d90e718"},
                "deepip": {"source_commit": "988564875d002b6d5d48d8dfb228cba3492dd776"},
                "models": {}, "environment": {},
            }, "mouse")

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
