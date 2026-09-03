import csv
import shlex
import tempfile
import unittest
from pathlib import Path

from rnaends2tracks.config import RunPlan
from rnaends2tracks.report import (
    _apa_b_gene_events, _browser_assets, _contrast_summary, _exact_funnel_rows,
    _html_table, _star_qc_rows, _top_apa_events, _top_dge_events,
    _top_enrichment_terms, _track_collections, _validate_ucsc_track_lines,
)


def write_tsv(path, fields, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


class ScientificReportTests(unittest.TestCase):
    def test_top_event_summaries_cover_dge_and_all_three_apa_methods(self):
        with tempfile.TemporaryDirectory() as temporary:
            results = Path(temporary)
            dge = results / "05_gene_expression" / "GRCm39" / "C4_primary_deseq2"
            dge_result = dge / "x.deseq2.tsv"
            write_tsv(dge_result, ["gene_id", "baseMean", "log2FoldChange", "pvalue", "padj"], [
                {"gene_id": "g1", "baseMean": "100", "log2FoldChange": "2", "pvalue": "0.001", "padj": "0.01"},
            ])
            write_tsv(dge / "result_index.tsv", ["contrast_id", "result_file"], [
                {"contrast_id": "x", "result_file": str(dge_result)},
            ])
            for method_dir, method, site_field, primary in (
                ("06_apa_a_mcell2019", "APA-A", "significant_sites", ""),
                ("06b_apa_a2_corrected", "APA-A2", "primary_sites", "true"),
                ("07_apa_b", "APA-B", "confirmed_sites", ""),
            ):
                root = results / method_dir / "GRCm39"
                summary = root / "stats" / "x.gene.tsv"
                write_tsv(summary, ["gene_id", "gene_padj", "shift", "max_abs_delta_PAU", site_field, "primary_gene"], [
                    {"gene_id": "g1", "gene_padj": "0.02", "shift": "distal",
                     "max_abs_delta_PAU": "0.3", site_field: "1", "primary_gene": primary},
                ])
                write_tsv(root / "stats" / "result_index.tsv", ["contrast_id", "gene_summary_file"], [
                    {"contrast_id": "x", "gene_summary_file": str(summary)},
                ])
            self.assertEqual(_top_dge_events(results, 0.05)[0]["direction"], "up")
            self.assertEqual(
                {row["method"] for row in _top_apa_events(results, 0.05)},
                {"APA-A", "APA-A2", "APA-B"},
            )

    def test_validated_apa_b_events_use_the_drimseq_primary_index(self):
        with tempfile.TemporaryDirectory() as temporary:
            results = Path(temporary)
            root = results / "07_apa_b" / "GRCm39" / "drimseq"
            summary = root / "x.gene_apa_summary.tsv"
            write_tsv(summary, ["gene_id", "gene_padj", "shift", "confirmed_sites"], [
                {"gene_id": "g1", "gene_padj": "0.01", "shift": "distal", "confirmed_sites": "2"},
            ])
            write_tsv(root / "result_index.tsv", ["contrast_id", "gene_summary_file"], [
                {"contrast_id": "x", "gene_summary_file": str(summary)},
            ])
            self.assertEqual(_apa_b_gene_events(results, 0.05)[0]["gene_id"], "g1")

    def test_top_enrichment_retains_ora_and_gsea_terms_independently(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            ora = root / "ora.tsv"; gsea = root / "gsea.tsv"
            fields = ["query", "database", "term_id", "term_name", "overlap_count", "NES", "padj"]
            write_tsv(ora, fields, [
                {"term_id": "ora1", "term_name": "ORA", "overlap_count": "5", "padj": "0.001"},
            ])
            write_tsv(gsea, fields, [
                {"term_id": "gsea1", "term_name": "GSEA", "NES": "2.1", "padj": "0.002"},
            ])
            terms = _top_enrichment_terms([{
                "analysis_type": "dge", "genome": "GRCm39", "contrast_id": "x",
                "ora_file": str(ora), "gsea_file": str(gsea),
            }], 0.05, limit_per_job=1)
            self.assertEqual({row["method"] for row in terms}, {"ORA", "GSEA"})

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
            gene_summary = apa_dir / f"{contrast_id}.gene_apa_summary.tsv"
            write_tsv(apa_result, ["pas_id", "padj"], [
                {"pas_id": "p1", "padj": "0.01"},
                {"pas_id": "p2", "padj": "0.2"},
            ])
            write_tsv(shift_result, ["gene_id", "shift"], [
                {"gene_id": "g1", "shift": "distal"},
                {"gene_id": "g2", "shift": "proximal"},
                {"gene_id": "g3", "shift": "no_shift"},
            ])
            write_tsv(gene_summary, ["gene_id", "gene_padj"], [
                {"gene_id": "g1", "gene_padj": "0.01"},
                {"gene_id": "g2", "gene_padj": "0.2"},
            ])
            write_tsv(
                apa_dir / "result_index.tsv",
                ["contrast_id", "result_file", "shift_file", "gene_summary_file", "tested_sites", "significant_sites"],
                [{"contrast_id": contrast_id, "result_file": str(apa_result),
                  "shift_file": str(shift_result), "gene_summary_file": str(gene_summary),
                  "tested_sites": "2", "significant_sites": "1"}],
            )
            write_tsv(
                results / "06_apa_a_mcell2019" / "GRCm39" / "candidate_pcpa.tsv",
                ["contrast_id", "pas_id"], [{"contrast_id": contrast_id, "pas_id": "p1"}],
            )

            apa_a2_dir = results / "06b_apa_a2_corrected" / "GRCm39" / "dexseq_a2"
            apa_a2_sites = apa_a2_dir / f"{contrast_id}.apa_a2_sites.tsv"
            apa_a2_genes = apa_a2_dir / f"{contrast_id}.apa_a2_genes.tsv"
            write_tsv(apa_a2_sites, ["pas_id", "padj", "primary_site"], [
                {"pas_id": "p1", "padj": "0.01", "primary_site": "true"},
                {"pas_id": "p2", "padj": "0.02", "primary_site": "false"},
                {"pas_id": "p3", "padj": "0.5", "primary_site": "false"},
            ])
            write_tsv(apa_a2_genes, ["gene_id", "gene_padj", "shift", "primary_gene"], [
                {"gene_id": "g1", "gene_padj": "0.01", "shift": "distal", "primary_gene": "true"},
                {"gene_id": "g2", "gene_padj": "0.02", "shift": "no_shift", "primary_gene": "false"},
            ])
            write_tsv(
                apa_a2_dir / "result_index.tsv",
                ["contrast_id", "result_file", "gene_summary_file", "tested_sites", "significant_sites", "primary_sites"],
                [{"contrast_id": contrast_id, "result_file": str(apa_a2_sites),
                  "gene_summary_file": str(apa_a2_genes), "tested_sites": "3",
                  "significant_sites": "2", "primary_sites": "1"}],
            )
            write_tsv(
                results / "06b_apa_a2_corrected" / "GRCm39" / "candidate_pcpa.tsv",
                ["contrast_id", "pas_id"], [{"contrast_id": contrast_id, "pas_id": "p1"}],
            )

            row = _contrast_summary(plan, results, [contrast])[0]
            self.assertEqual(row["dge_tested_genes"], 3)
            self.assertEqual((row["dge_significant"], row["dge_up"], row["dge_down"]), (2, 1, 1))
            self.assertEqual((row["apa_a_tested_sites"], row["apa_a_significant_sites"]), (2, 1))
            self.assertEqual((row["apa_a_tested_genes"], row["apa_a_significant_genes"]), (2, 1))
            self.assertEqual((row["apa_a_distal_genes"], row["apa_a_proximal_genes"]), (1, 1))
            self.assertEqual(row["apa_a_pcpa"], 1)
            self.assertEqual((row["apa_a2_tested_sites"], row["apa_a2_significant_sites"]), (3, 2))
            self.assertEqual((row["apa_a2_primary_sites"], row["apa_a2_primary_genes"]), (1, 1))
            self.assertEqual((row["apa_a2_distal_genes"], row["apa_a2_proximal_genes"]), (1, 0))
            self.assertEqual(row["apa_a2_pcpa"], 1)

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

    def test_browser_assets_are_grouped_and_use_one_line_ucsc_syntax(self):
        with tempfile.TemporaryDirectory() as temporary:
            results = Path(temporary) / "results"
            outdir = results / "10_reports"
            outdir.mkdir(parents=True)
            for group in ("all_reads/raw", "active_pas/robust_cpm"):
                folder = results / "09_tracks" / group
                folder.mkdir(parents=True)
                for strand in ("plus", "minus"):
                    token = group.replace("/", ".")
                    (folder / f"sample.{token}.transcript_{strand}.bw").write_bytes(b"bw")
            project = {"tracks": {
                "ucsc_bigdata_url_prefix": "http://example.test/tracks",
                "ucsc_negate_minus_tracks": True,
                "ucsc_view_limits": "0:12",
            }}
            plan = RunPlan(project, [], [], [], {}, {})

            outputs, collections = _browser_assets(plan, results, outdir)

            self.assertEqual([row["collection"] for row in collections], [
                "all_reads/raw", "active_pas/robust_cpm",
            ])
            inventory = (outdir / "bigwig_collections.txt").read_text(encoding="utf-8")
            self.assertIn("[all_reads/raw]", inventory)
            self.assertNotIn("\t", inventory)
            combined = (
                outdir / "ucsc_track_descriptors" / "UCSC_bigWig_tracks.oneline.txt"
            ).read_text(encoding="utf-8").splitlines()
            track_lines = [line for line in combined if line.startswith("track ")]
            self.assertEqual(len(track_lines), 4)
            self.assertTrue(all(" type=bigWig " in line for line in track_lines))
            self.assertTrue(all("bigDataUrl=http://example.test/tracks/sample." in line for line in track_lines))
            self.assertTrue(all("/all_reads/" not in line and "/active_pas/" not in line for line in track_lines))
            self.assertEqual(sum("negateValues=on" in line for line in track_lines), 2)
            self.assertTrue(all("viewLimits=0:12" in line for line in track_lines))
            descriptor_dir = outdir / "ucsc_track_descriptors"
            self.assertIn(descriptor_dir / "all_reads.txt", outputs)
            self.assertIn(descriptor_dir / "active_pas.txt", outputs)
            self.assertIn(descriptor_dir / "UCSC_trackDb.txt", outputs)
            self.assertIn(descriptor_dir / "UCSC_descriptor_validation.tsv", outputs)
            self.assertTrue(all(len(line) == len(line.rstrip("\r\n")) for line in track_lines))
            self.assertTrue(all(
                len(next(token for token in shlex.split(line)
                         if token.startswith("description=")).split("=", 1)[1]) <= 60
                for line in track_lines
            ))

    def test_ucsc_validator_rejects_markdown_url_and_multiline_records(self):
        valid = (
            'track type=bigWig name="rna_ends_0001" description="valid" '
            'bigDataUrl=http://example.test/sample.bw visibility=full color=0,102,204 '
            'viewLimits=0:12'
        )
        self.assertEqual(_validate_ucsc_track_lines([valid]), 1)
        with self.assertRaisesRegex(RuntimeError, "bigDataUrl"):
            _validate_ucsc_track_lines([valid.replace(
                "http://example.test/sample.bw", "[http://example.test](http://example.test)/sample.bw",
            )])
        with self.assertRaisesRegex(RuntimeError, "not one line"):
            _validate_ucsc_track_lines([valid + "\ntrack type=bigWig"])

    def test_track_collection_counts_strands(self):
        with tempfile.TemporaryDirectory() as temporary:
            results = Path(temporary)
            folder = results / "09_tracks" / "filtered_ends" / "cpm"
            folder.mkdir(parents=True)
            (folder / "s1.transcript_plus.bw").write_bytes(b"bw")
            (folder / "s1.transcript_minus.bw").write_bytes(b"bw")
            _, rows = _track_collections(results)
            self.assertEqual(rows[0]["bigwigs"], 2)
            self.assertEqual(rows[0]["transcript_plus"], 1)
            self.assertEqual(rows[0]["transcript_minus"], 1)

    def test_qc_sources_are_summarized_without_reinterpreting_counts(self):
        with tempfile.TemporaryDirectory() as temporary:
            results = Path(temporary)
            star = results / "02_alignment" / "S1" / "lanes" / "S1.T01.L001.star.Log.final.out"
            star.parent.mkdir(parents=True)
            star.write_text(
                "Number of input reads | 100\nUniquely mapped reads number | 80\n"
                "Uniquely mapped reads % | 80.00%\n% of reads mapped to multiple loci | 10.00%\n"
                "% of reads unmapped: too short | 5.00%\n",
                encoding="utf-8",
            )
            audit = results / "03_exact_ends" / "GRCm39" / "S1" / "end_audit.json"
            audit.parent.mkdir(parents=True)
            audit.write_text(
                '{"sample_id":"S1","genome":"GRCm39","C0":80,"C1":72,'
                '"C1S":8,"C2":60,"C2R":12,"mask_rescued":3}\n',
                encoding="utf-8",
            )
            star_rows = _star_qc_rows(results)
            funnel_rows = _exact_funnel_rows(results)
            self.assertEqual(star_rows[0]["uniquely_mapped_reads"], "80")
            self.assertEqual(star_rows[0]["uniquely_mapped_pct"], "80.00%")
            self.assertEqual(funnel_rows[0]["C0"], 80)
            self.assertEqual(funnel_rows[0]["C1_over_C0_pct"], "90.00")
            self.assertEqual(funnel_rows[0]["C2_over_C1_pct"], "83.33")


if __name__ == "__main__":
    unittest.main()
