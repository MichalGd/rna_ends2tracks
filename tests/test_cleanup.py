import json
import tempfile
import unittest
from pathlib import Path

from rnaends2tracks import __version__
from rnaends2tracks.cleanup import clean_intermediates
from rnaends2tracks.config import RunPlan
from rnaends2tracks.receipts import sha256

REQUIRED = ["02_alignment", "03_exact_ends", "04_active_pas", "05_gene_expression",
            "06_apa_a_mcell2019", "09_tracks", "10_reports"]


def _plan() -> RunPlan:
    return RunPlan(
        project={
            "cleanup": {"enabled": True, "keep_trimmed_fastq": False, "keep_lane_bams": False,
                "keep_apa_sample_extraction": False, "keep_track_strand_bams": False,
                "keep_track_bedgraphs": False},
            "modules": {"gene_expression": True, "apa_a": True, "tracks": True},
            "apa_b": {"enabled": False},
        },
        samples=[{"sample_id": "S1"}],
        sample_rows=[{"sample_id": "S1", "technical_replicate_id": "T01", "lane_id": "L001"}],
        contrasts=[], reference={},
    )


def _write_success_receipts(results: Path, workflow_version: str = __version__) -> None:
    for module in REQUIRED:
        module_dir = results / module
        output = module_dir / "final.output"
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_bytes(b"complete")
        (module_dir / "run_receipt.json").write_text(json.dumps({
            "workflow_version": workflow_version, "exit_status": 0,
            "outputs": [{
                "path": str(output.resolve()), "size": output.stat().st_size,
                "sha256": sha256(output),
            }],
        }), encoding="utf-8")


class CleanupTests(unittest.TestCase):
    def test_default_cleanup_removes_only_allowlisted_intermediates(self):
        with tempfile.TemporaryDirectory() as temporary:
            results = Path(temporary) / "project" / "results"
            trimmed = results / "01_qc" / "trimmed_fastq" / "S1.trimmed.fastq.gz"
            unrelated = results / "01_qc" / "trimmed_fastq" / "README.txt"
            lane = results / "02_alignment" / "S1" / "lanes" / "S1.T01.L001.bam"
            star = results / "02_alignment" / "S1" / "lanes" / "S1.T01.L001.star.Aligned.out.bam"
            manual = results / "02_alignment" / "S1" / "lanes" / "manual.bam"
            final_bam = results / "02_alignment" / "S1" / "S1.bam"
            track_dir = results / "09_tracks" / ".intermediate" / "S1"
            strand_bam = track_dir / "S1.all_reads.raw.plus.strand.bam"
            bedgraph = track_dir / "S1.C2.cpm.plus.bedGraph"
            final_bigwig = results / "09_tracks" / "filtered_ends" / "cpm" / "S1.bw"
            final_matrix = results / "04_active_pas" / "GRCm39" / "C3_active_pas_counts.tsv"
            for path in (trimmed, unrelated, lane, star, manual, final_bam, strand_bam,
                         bedgraph, final_bigwig, final_matrix):
                path.parent.mkdir(parents=True, exist_ok=True); path.write_bytes(b"test")
            _write_success_receipts(results)

            manifest = clean_intermediates(_plan(), results)

            self.assertIsNotNone(manifest)
            self.assertFalse(trimmed.exists()); self.assertFalse(lane.exists()); self.assertFalse(star.exists())
            self.assertFalse(strand_bam.exists()); self.assertFalse(bedgraph.exists())
            self.assertTrue(unrelated.is_file()); self.assertTrue(manual.is_file()); self.assertTrue(final_bam.is_file())
            self.assertTrue(final_bigwig.is_file()); self.assertTrue(final_matrix.is_file())
            text = manifest.read_text(encoding="utf-8")
            for category in ("trimmed_fastq", "lane_alignment_bam", "track_strand_bam", "track_bedgraph"):
                self.assertIn(category, text)

    def test_cleanup_refuses_incomplete_workflow(self):
        with tempfile.TemporaryDirectory() as temporary:
            results = Path(temporary) / "project" / "results"
            trimmed = results / "01_qc" / "trimmed_fastq" / "S1.trimmed.fastq.gz"
            trimmed.parent.mkdir(parents=True, exist_ok=True); trimmed.write_bytes(b"test")
            with self.assertRaisesRegex(RuntimeError, "successful complete workflow"):
                clean_intermediates(_plan(), results)
            self.assertTrue(trimmed.is_file())

    def test_hotfix_cleanup_accepts_audited_alpha9_receipts(self):
        with tempfile.TemporaryDirectory() as temporary:
            results = Path(temporary) / "project" / "results"
            trimmed = results / "01_qc" / "trimmed_fastq" / "S1.trimmed.fastq.gz"
            trimmed.parent.mkdir(parents=True, exist_ok=True); trimmed.write_bytes(b"test")
            _write_success_receipts(results, workflow_version="0.1.0a9")
            clean_intermediates(_plan(), results)
            self.assertFalse(trimmed.exists())

    def test_repeated_cleanup_preserves_manifest(self):
        with tempfile.TemporaryDirectory() as temporary:
            results = Path(temporary) / "project" / "results"
            trimmed = results / "01_qc" / "trimmed_fastq" / "S1.trimmed.fastq.gz"
            trimmed.parent.mkdir(parents=True, exist_ok=True); trimmed.write_bytes(b"test")
            _write_success_receipts(results)
            manifest = clean_intermediates(_plan(), results)
            first = manifest.read_text(encoding="utf-8")
            clean_intermediates(_plan(), results, force=True)
            self.assertEqual(manifest.read_text(encoding="utf-8"), first)


if __name__ == "__main__":
    unittest.main()
