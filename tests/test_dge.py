import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from rnaends2tracks.config import RunPlan
from rnaends2tracks.dge import gene_expression


class GeneExpressionTests(unittest.TestCase):
    def test_genome_output_directory_exists_before_featurecounts(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            results = root / "results"
            reference = {"assembly": "GRCm39", "gtf": str(root / "genes.gtf")}
            plan = RunPlan(
                project={
                    "project_id": "dge_test",
                    "modules": {"gene_expression": True},
                    "design": "~ condition",
                    "reporting": {"fdr": 0.05},
                    "resources": {"dge": {
                        "featurecounts_threads": 1,
                        "contrast_parallel_jobs": 1,
                        "contrast_memory_gb": 1,
                    }},
                },
                samples=[{"sample_id": "S1", "genome": "GRCm39"}],
                sample_rows=[], contrasts=[], reference=reference,
                references={"GRCm39": reference},
            )

            def fake_run(command, _log, dry_run=False, cwd=None, env=None):
                self.assertFalse(dry_run)
                if command[0] == "featureCounts":
                    target = Path(command[command.index("-o") + 1])
                    self.assertTrue(target.parent.is_dir())

            with (
                patch("rnaends2tracks.dge.require_tools"),
                patch("rnaends2tracks.dge.signature_for", return_value="signature"),
                patch("rnaends2tracks.dge.receipt_valid", return_value=False),
                patch("rnaends2tracks.dge.run", side_effect=fake_run),
                patch("rnaends2tracks.dge.write_receipt"),
                patch("rnaends2tracks.dge._write_c4_c5_diagnostic"),
            ):
                gene_expression(plan, results, root / "scripts")


if __name__ == "__main__":
    unittest.main()
