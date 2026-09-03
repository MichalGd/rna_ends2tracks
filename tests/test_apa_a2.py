import csv
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from rnaends2tracks.apa_a2 import apa_a2
from rnaends2tracks.config import RunPlan


class ApaA2StageTests(unittest.TestCase):
    def test_stage_has_independent_outputs_receipt_and_r_script(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            results = root / "results"
            active = results / "04_active_pas" / "GRCm39"
            active.mkdir(parents=True)
            counts = active / "C3_active_pas_counts.tsv"
            catalog = active / "active_pas_catalog.tsv"
            counts.write_text("pas_id\n", encoding="utf-8")
            catalog.write_text("pas_id\n", encoding="utf-8")
            contrast = {
                "contrast_id": "treated_vs_control",
                "genome": "GRCm39",
                "denominator": "control",
                "numerator": "treated",
            }
            samples = [
                {"sample_id": "control_1", "genome": "GRCm39"},
                {"sample_id": "treated_1", "genome": "GRCm39"},
            ]
            project = {
                "design": "~ condition",
                "reporting": {"fdr": 0.05, "min_abs_delta_pau": 0.10},
                "resources": {"apa_a2": {"contrast_parallel_jobs": 2, "contrast_memory_gb": 16}},
            }
            reference = {"assembly": "GRCm39"}
            plan = RunPlan(project, samples, samples, [contrast], reference, {"GRCm39": reference})
            captured = {}

            def fake_run_r_contrasts(**kwargs):
                captured.update(kwargs)
                outdir = kwargs["outdir"]
                outdir.mkdir(parents=True, exist_ok=True)
                contrast_id = contrast["contrast_id"]
                paths = {
                    "result_file": outdir / f"{contrast_id}.apa_a2_sites.tsv",
                    "gene_summary_file": outdir / f"{contrast_id}.apa_a2_genes.tsv",
                    "pair_delta_file": outdir / f"{contrast_id}.apa_a2_pair_deltas.tsv",
                    "shift_file": outdir / f"{contrast_id}.apa_a2_shifts.tsv",
                    "audit_file": outdir / f"{contrast_id}.apa_a2_audit.tsv",
                }
                for path in paths.values():
                    path.write_text("field\n", encoding="utf-8")
                with kwargs["index_path"].open("w", encoding="utf-8", newline="") as handle:
                    writer = csv.DictWriter(
                        handle,
                        fieldnames=["contrast_id", *paths],
                        delimiter="\t",
                        lineterminator="\n",
                    )
                    writer.writeheader()
                    writer.writerow({"contrast_id": contrast_id, **paths})

            def fake_pcpa(_candidate, _index, output, _fdr, _min_delta):
                output.write_text("pas_id\tcontrast_id\n", encoding="utf-8")

            with (
                patch("rnaends2tracks.apa_a2.require_tools"),
                patch("rnaends2tracks.apa_a2.run_r_contrasts", side_effect=fake_run_r_contrasts),
                patch("rnaends2tracks.apa_a2._write_differential_pcpa", side_effect=fake_pcpa),
            ):
                apa_a2(plan, results, root / "scripts")

            self.assertEqual(captured["script"].name, "dexseq_all_pairs_a2.R")
            self.assertEqual(captured["module"], "apa_a2_GRCm39")
            self.assertIn(".apa_a2_audit.tsv", captured["output_suffixes"])
            receipt = json.loads(
                (results / "06b_apa_a2_corrected" / "run_receipt.json").read_text(encoding="utf-8")
            )
            self.assertEqual(receipt["module"], "apa_a2_corrected")
            self.assertFalse((results / "06_apa_a_mcell2019").exists())


if __name__ == "__main__":
    unittest.main()
