import csv
import tempfile
import unittest
from pathlib import Path

from rnaends2tracks.config import RunPlan
from rnaends2tracks.report import _contrast_summary, _html_table


def write_tsv(path, fields, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


class ScientificReportTests(unittest.TestCase):
    def test_contrast_summary_recounts_dge_apa_and_pcpa_sources(self):
        with tempfile.TemporaryDirectory() as temporary:
            results = Path(temporary)
            contrast_id = "treated_vs_control"
            contrast = {
                "genome": "GRCm39", "contrast_id": contrast_id,
                "numerator": "treated", "denominator": "control", "design_mode": "paired",
            }
            reference = {"assembly": "GRCm39"}
            plan = RunPlan(
                {"reporting": {"fdr": 0.05}}, [], [], [contrast], reference,
                {"GRCm39": reference},
            )

            dge_dir = results / "05_gene_expression" / "GRCm39" / "C4_primary_deseq2"
            dge_result = dge_dir / f"{contrast_id}.deseq2.tsv"
            write_tsv(dge_result, ["gene_id", "log2FoldChange", "padj"], [
                {"gene_id": "g1", "log2FoldChange": "2", "padj": "0.01"},
                {"gene_id": "g2", "log2FoldChange": "-1", "padj": "0.02"},
                {"gene_id": "g3", "log2FoldChange": "3", "padj": "0.5"},
            ])
            write_tsv(dge_dir / "result_index.tsv", ["contrast_id", "result_file", "significant"], [{
                "contrast_id": contrast_id, "result_file": str(dge_result), "significant": "2",
            }])

            apa_dir = results / "06_apa_a_mcell2019" / "GRCm39" / "dexseq"
            apa_result = apa_dir / f"{contrast_id}.dexseq.tsv"
            shift_result = apa_dir / f"{contrast_id}.apa_shift.tsv"
            write_tsv(apa_result, ["pas_id", "padj"], [
                {"pas_id": "p1", "padj": "0.01"},
                {"pas_id": "p2", "padj": "0.2"},
            ])
            write_tsv(shift_result, ["gene_id", "shift"], [
                {"gene_id": "g1", "shift": "distal"},
                {"gene_id": "g2", "shift": "proximal"},
                {"gene_id": "g3", "shift": "no_shift"},
            ])
            write_tsv(
                apa_dir / "result_index.tsv",
                ["contrast_id", "result_file", "shift_file", "tested_sites", "significant_sites"],
                [{"contrast_id": contrast_id, "result_file": str(apa_result),
                  "shift_file": str(shift_result), "tested_sites": "2", "significant_sites": "1"}],
            )
            write_tsv(
                results / "06_apa_a_mcell2019" / "GRCm39" / "candidate_pcpa.tsv",
                ["contrast_id", "pas_id"], [{"contrast_id": contrast_id, "pas_id": "p1"}],
            )

            row = _contrast_summary(plan, results, [contrast])[0]
            self.assertEqual(row["dge_tested_genes"], 3)
            self.assertEqual((row["dge_significant"], row["dge_up"], row["dge_down"]), (2, 1, 1))
            self.assertEqual((row["apa_a_tested_sites"], row["apa_a_significant_sites"]), (2, 1))
            self.assertEqual((row["apa_a_distal_genes"], row["apa_a_proximal_genes"]), (1, 1))
            self.assertEqual(row["apa_a_pcpa"], 1)

    def test_contrast_summary_rejects_stale_index_counts(self):
        with tempfile.TemporaryDirectory() as temporary:
            results = Path(temporary)
            contrast = {"genome": "GRCm39", "contrast_id": "a_vs_b"}
            reference = {"assembly": "GRCm39"}
            plan = RunPlan(
                {"reporting": {"fdr": 0.05}}, [], [], [contrast], reference,
                {"GRCm39": reference},
            )
            dge_dir = results / "05_gene_expression" / "GRCm39" / "C4_primary_deseq2"
            result = dge_dir / "a_vs_b.deseq2.tsv"
            write_tsv(result, ["gene_id", "log2FoldChange", "padj"], [
                {"gene_id": "g1", "log2FoldChange": "1", "padj": "0.01"},
            ])
            write_tsv(dge_dir / "result_index.tsv", ["contrast_id", "result_file", "significant"], [{
                "contrast_id": "a_vs_b", "result_file": str(result), "significant": "0",
            }])

            with self.assertRaisesRegex(RuntimeError, "DGE summary mismatch"):
                _contrast_summary(plan, results, [contrast])

    def test_html_table_escapes_source_values(self):
        rendered = _html_table([{"name": "<unsafe>"}], ["name"])
        self.assertIn("&lt;unsafe&gt;", rendered)
        self.assertNotIn("<unsafe>", rendered)


if __name__ == "__main__":
    unittest.main()
