import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from rnaends2tracks.config import RunPlan
from rnaends2tracks.tracks import (
    _sample_tracks_subset,
    make_c0_tracks,
    make_c0_tracks_for_sample,
    make_tracks,
)


class EarlyTrackTests(unittest.TestCase):
    def test_sample_ready_c0_tracks_use_only_raw_and_cpm(self):
        sample = {"sample_id": "S1", "genome": "GRCm39"}
        reference = {"assembly": "GRCm39", "chrom_sizes": "chrom.sizes"}
        project = {
            "modules": {"tracks": True},
            "tracks": {
                "early_c0": True,
                "families": {"all_reads": True},
                "normalizations": {
                    "raw": True, "cpm": True, "deseq2": True, "robust_cpm": True,
                },
                "generate_bigwigs": False,
            },
        }
        plan = RunPlan(project, [sample], [], [], reference, {"GRCm39": reference})

        with (
            patch("rnaends2tracks.tracks.require_tools"),
            patch("rnaends2tracks.tracks._sample_tracks_subset", return_value=([], [])) as subset,
        ):
            make_c0_tracks_for_sample(plan, Path("results"), sample)

        self.assertEqual(subset.call_args.args[4], ("all_reads",))
        self.assertEqual(subset.call_args.args[5], ("raw", "cpm"))
        self.assertEqual(subset.call_args.args[6], "tracks_c0")

    def test_final_track_stage_reuses_early_c0_receipt(self):
        with tempfile.TemporaryDirectory() as temporary:
            results = Path(temporary) / "results"
            bam = results / "02_alignment" / "S1" / "S1.bam"
            bam.parent.mkdir(parents=True)
            bam.write_bytes(b"fixture")
            sizes = Path(temporary) / "chrom.sizes"
            sizes.write_text("chr1\t100\n", encoding="utf-8")
            sample = {"sample_id": "S1", "genome": "GRCm39"}
            reference = {"assembly": "GRCm39", "chrom_sizes": str(sizes)}
            project = {
                "modules": {"tracks": True},
                "tracks": {
                    "early_c0": True,
                    "families": {
                        "all_reads": True, "exact_ends": False, "filtered_ends": False,
                        "rejected_ends": False, "active_pas": False,
                    },
                    "normalizations": {
                        "raw": True, "cpm": False, "deseq2": False, "robust_cpm": False,
                    },
                    "generate_bigwigs": True,
                    "retain_bedgraph": False,
                },
                "resources": {"tracks": {"parallel_jobs": 1, "samtools_threads": 2, "memory_gb": 4}},
            }
            plan = RunPlan(project, [sample], [], [], reference, {"GRCm39": reference})
            deliverable = results / "09_tracks" / "all_reads" / "raw" / "S1.raw.bw"
            row = {
                "sample_id": "S1", "genome": "GRCm39", "family": "all_reads",
                "normalization": "raw", "denominator": 100, "scale": "1", "count_universe": "C0",
            }

            def fake_subset(*_args, **_kwargs):
                deliverable.parent.mkdir(parents=True, exist_ok=True)
                deliverable.write_bytes(b"bigwig")
                return [deliverable], [row]

            with (
                patch("rnaends2tracks.tracks.require_tools"),
                patch("rnaends2tracks.tracks._run_track_subset", side_effect=fake_subset) as subset_mock,
            ):
                make_c0_tracks(plan, results)
                make_tracks(plan, results)

            self.assertEqual(subset_mock.call_count, 1)
            self.assertTrue((results / "09_tracks" / "run_receipt.json").is_file())
            combined = (results / "09_tracks" / "track_normalization.tsv").read_text()
            self.assertIn("S1\tGRCm39\tall_reads\traw", combined)

    def test_c0_raw_and_cpm_reuse_one_bam_per_strand(self):
        with tempfile.TemporaryDirectory() as temporary:
            results = Path(temporary) / "results"
            bam = results / "02_alignment" / "S1" / "S1.bam"
            bam.parent.mkdir(parents=True)
            bam.write_bytes(b"fixture")
            sizes = Path(temporary) / "chrom.sizes"
            sizes.write_text("chr1\t100\n", encoding="utf-8")
            project = {
                "tracks": {
                    "early_c0": True,
                    "families": {
                        "all_reads": True, "exact_ends": False, "filtered_ends": False,
                        "rejected_ends": False, "active_pas": False,
                    },
                    "normalizations": {
                        "raw": True, "cpm": True, "deseq2": False, "robust_cpm": False,
                    },
                    "generate_bigwigs": False,
                    "retain_bedgraph": True,
                },
                "resources": {"tracks": {"parallel_jobs": 1, "samtools_threads": 2, "memory_gb": 4}},
            }
            sample = {"sample_id": "S1", "genome": "GRCm39"}
            reference = {"assembly": "GRCm39", "chrom_sizes": str(sizes)}
            plan = RunPlan(project, [sample], [], [], reference, {"GRCm39": reference})

            def fake_strand_bam(_bam, _strand, output, _threads, _log):
                output.parent.mkdir(parents=True, exist_ok=True)
                output.write_bytes(b"strand")

            def fake_coverage(_bam, _strand, output, _scale, _log):
                output.parent.mkdir(parents=True, exist_ok=True)
                output.write_text("chr1\t0\t1\t1\n", encoding="utf-8")

            with (
                patch("rnaends2tracks.tracks._bam_count", return_value=100),
                patch("rnaends2tracks.tracks._strand_bam", side_effect=fake_strand_bam) as strand_mock,
                patch("rnaends2tracks.tracks._all_read_bedgraph", side_effect=fake_coverage) as coverage_mock,
            ):
                outputs, rows = _sample_tracks_subset(
                    plan, results, sample, False, ("all_reads",), ("raw", "cpm"), "tracks_c0"
                )

            self.assertEqual(strand_mock.call_count, 2)
            self.assertEqual(coverage_mock.call_count, 4)
            self.assertEqual({row["normalization"] for row in rows}, {"raw", "cpm"})
            self.assertEqual(len([path for path in outputs if path.suffix == ".bedGraph"]), 4)


if __name__ == "__main__":
    unittest.main()
