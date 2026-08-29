import csv
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from rnaends2tracks.config import RunPlan
from rnaends2tracks.enrichment import (
    _apa_gene_table, _dge_gene_table, _ensure_apa_a_gene_summaries,
)
from rnaends2tracks.provenance import generate_provenance_dashboard


def write_tsv(path, fields, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)


def read_tsv(path):
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


class EnrichmentTests(unittest.TestCase):
    def test_legacy_apa_a_index_is_regenerated_exactly_before_enrichment(self):
        with tempfile.TemporaryDirectory() as temporary:
            results = Path(temporary) / "results"
            index = results / "06_apa_a_mcell2019" / "GRCm39" / "dexseq" / "result_index.tsv"
            write_tsv(index, ["contrast_id", "result_file"], [{"contrast_id": "x", "result_file": "old.tsv"}])
            plan = RunPlan({}, [], [], [], {}, {"GRCm39": {"assembly": "GRCm39"}})

            def regenerate(*_args, **_kwargs):
                summary = index.parent / "x.gene_apa_summary.tsv"
                write_tsv(summary, ["gene_id", "gene_padj"], [{"gene_id": "g1", "gene_padj": "0.1"}])
                write_tsv(index, ["contrast_id", "gene_summary_file"], [
                    {"contrast_id": "x", "gene_summary_file": str(summary)},
                ])

            with patch("rnaends2tracks.apa_mcell.apa_statistics_stage", side_effect=regenerate) as mocked:
                _ensure_apa_a_gene_summaries(plan, results, Path(temporary))
            mocked.assert_called_once()

    def test_dge_preparation_uses_tested_background_and_directional_foregrounds(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary); source = root / "dge.tsv"; target = root / "prepared.tsv"
            write_tsv(source, ["gene_id", "log2FoldChange", "stat", "pvalue", "padj"], [
                {"gene_id": "g1", "log2FoldChange": "2", "stat": "5", "pvalue": "1e-6", "padj": "0.01"},
                {"gene_id": "g2", "log2FoldChange": "-3", "stat": "-4", "pvalue": "1e-5", "padj": "0.02"},
                {"gene_id": "g3", "log2FoldChange": "0.3", "stat": "1", "pvalue": "0.2", "padj": "0.5"},
            ])
            self.assertEqual(_dge_gene_table(source, target, 0.05, 1.0), (3, 2))
            rows = read_tsv(target)
            self.assertEqual(rows[0]["direction"], "any_dge;upregulated")
            self.assertEqual(rows[1]["direction"], "any_dge;downregulated")
            self.assertEqual(rows[2]["foreground"], "0")

    def test_apa_preparation_separates_shift_and_pcpa_queries(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary); source = root / "apa.tsv"; pcpa = root / "pcpa.tsv"; target = root / "prepared.tsv"
            write_tsv(source, ["gene_id", "gene_padj", "max_abs_delta_PAU", "signed_shift_score", "shift"], [
                {"gene_id": "g1", "gene_padj": "0.01", "max_abs_delta_PAU": "0.3", "signed_shift_score": "0.3", "shift": "distal"},
                {"gene_id": "g2", "gene_padj": "0.5", "max_abs_delta_PAU": "0.2", "signed_shift_score": "-0.2", "shift": "proximal"},
            ])
            write_tsv(pcpa, ["gene_id", "contrast_id", "delta_PAU"], [
                {"gene_id": "g2", "contrast_id": "t_vs_c", "delta_PAU": "-0.4"},
            ])
            self.assertEqual(_apa_gene_table(source, pcpa, target, "t_vs_c", 0.05, 0.1), (2, 2))
            rows = {row["gene_id"]: row for row in read_tsv(target)}
            self.assertEqual(rows["g1"]["direction"], "any_apa;distal_shift")
            self.assertEqual(rows["g2"]["direction"], "any_pcpa;pcpa_decreased")

    def test_provenance_dashboard_inventories_receipts_and_outputs(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary); results = root / "results"; outdir = results / "10_reports"
            (results / "module").mkdir(parents=True); outdir.mkdir(parents=True)
            (results / "module" / "value.tsv").write_text("a\n1\n", encoding="utf-8")
            (results / "module" / "run_receipt.json").write_text(json.dumps({
                "module": "test", "exit_status": 0, "workflow_version": "x", "outputs": [],
            }), encoding="utf-8")
            config = root / "config.conf"; sheet = root / "samples.csv"
            config.write_text("PROJECT_ID=x\n", encoding="utf-8"); sheet.write_text("sample_id\n", encoding="utf-8")
            reference = {"assembly": "GRCm39", "species": "mouse", "release": "vM31"}
            plan = RunPlan({"project_id": "x", "_config_path": str(config), "_samplesheet_path": str(sheet),
                            "apa_b": {"enabled": False}}, [], [], [], reference, {"GRCm39": reference})
            outputs = generate_provenance_dashboard(plan, results, outdir)
            self.assertTrue(all(path.is_file() for path in outputs))
            self.assertIn("module/value.tsv", (outdir / "provenance_dashboard" / "output_manifest.tsv").read_text())
            self.assertEqual(json.loads((outdir / "provenance_dashboard" / "dashboard.json").read_text())
                             ["apa_b_interpretation_status"], "DISABLED_NOT_VALIDATED")


if __name__ == "__main__":
    unittest.main()
