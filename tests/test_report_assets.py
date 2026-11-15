import base64
import json
import shutil
import tempfile
import unittest
from html.parser import HTMLParser
from pathlib import Path
from unittest.mock import patch
from urllib.parse import unquote

from rnaends2tracks.config import RunPlan
from rnaends2tracks.report import _image_gallery, make_report
from rnaends2tracks.report_assets import copy_report_image, repair_report_figures


PNG = base64.b64decode(
    "iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAQAAAC1HAwCAAAAC0lEQVR42mP8/x8AAwMCAO+aB9sAAAAASUVORK5CYII="
)
SVG = b'<svg xmlns="http://www.w3.org/2000/svg" width="200" height="100"><rect width="200" height="100" fill="blue"/></svg>'


class Images(HTMLParser):
    def __init__(self, contents):
        super().__init__()
        self.sources = []
        self.links = []
        self.feed(contents)

    def handle_starttag(self, tag, attrs):
        attributes = dict(attrs)
        if tag == "img":
            self.sources.append(attributes["src"])
        elif tag == "a":
            self.links.append(attributes.get("href", ""))


class ReportAssetTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.results = self.root / "results"
        self.outdir = self.results / "10_reports"
        self.outdir.mkdir(parents=True)
        self.pca = self.results / "05_gene_expression/GRCm39/C4_primary_deseq2/C4_vst_pca.png"
        self.pca.parent.mkdir(parents=True)
        self.pca.write_bytes(PNG)
        self.rseqc = self.results / "01_qc/rseqc/gene_body_coverage.svg"
        self.rseqc.parent.mkdir(parents=True)
        self.rseqc.write_bytes(SVG)

    def original(self, body):
        report = self.outdir / "report.html"
        report.write_bytes(("<!doctype html><html><body>" + body + "</body></html>").encode("utf-8"))
        return report

    def assert_downloadable(self, report):
        download = self.root / "download" / "10_reports"
        shutil.copytree(self.outdir, download)
        figures = Images((download / report.name).read_text(encoding="utf-8"))
        self.assertTrue(figures.sources)
        for value in figures.sources:
            path = (download / unquote(value)).resolve()
            self.assertTrue(path.is_relative_to(download.resolve()), value)
            self.assertTrue(path.is_file(), value)
        return figures

    def test_gallery_collects_png_and_svg_inside_downloaded_folder(self):
        records = []
        report = self.original(_image_gallery(self.outdir, [("PCA", self.pca), ("QC", self.rseqc)], records))
        figures = self.assert_downloadable(report)
        self.assertEqual(len(figures.sources), 2)
        self.assertEqual(set(figures.sources), set(figures.links))
        self.assertEqual(len(records), 2)

    def test_same_basename_different_contents_do_not_collide(self):
        other = self.results / "another/sample.png"
        first = self.results / "first/sample.png"
        other.parent.mkdir()
        first.parent.mkdir()
        other.write_bytes(PNG)
        first.write_bytes(PNG + b"different image")
        self.assertNotEqual(copy_report_image(self.outdir, other)["asset"],
                            copy_report_image(self.outdir, first)["asset"])

    def test_identical_images_share_an_asset_and_corrupted_copy_is_repaired(self):
        first = copy_report_image(self.outdir, self.pca)
        duplicate = self.pca.with_name("duplicate.png")
        duplicate.write_bytes(PNG)
        self.assertEqual(copy_report_image(self.outdir, duplicate)["asset"], first["asset"])
        target = self.outdir / first["asset"]
        target.write_bytes(b"corrupt")
        copy_report_image(self.outdir, self.pca)
        self.assertEqual(target.read_bytes(), PNG)

    def test_repair_rewrites_images_and_figure_links_but_preserves_tables_and_scripts(self):
        unchanged = '<table><tr><td>12 &amp; 34</td></tr></table>\r\n<script>const marker="<img src=bad>";</script>'
        report = self.original(
            unchanged + "<a href='../05_gene_expression/GRCm39/C4_primary_deseq2/C4_vst_pca.png'>"
            "<img alt='PCA &amp; samples' src='../05_gene_expression/GRCm39/C4_primary_deseq2/C4_vst_pca.png'></a>"
            "<img src='../01_qc/rseqc/gene_body_coverage.svg' alt='QC'>"
            "<a href='../05_gene_expression/counts.tsv'>Complete data</a>"
        )
        original_bytes = report.read_bytes()
        summary = repair_report_figures(report)
        output = Path(summary["html"])
        self.assertEqual(report.read_bytes(), original_bytes)
        contents = output.read_bytes().decode()
        self.assertIn(unchanged, contents)
        self.assertIn("href='../05_gene_expression/counts.tsv'", contents)
        figures = self.assert_downloadable(output)
        self.assertIn(figures.sources[0], figures.links)
        self.assertEqual(summary["local_images"], 2)

    def test_url_encoded_names_and_absolute_file_uris_are_supported(self):
        named = self.pca.with_name("space & figure.png")
        named.write_bytes(PNG)
        report = self.original(
            "<img src='../05_gene_expression/GRCm39/C4_primary_deseq2/space%20%26%20figure.png'/>"
            + '<img src="' + self.rseqc.as_uri() + '">'
        )
        summary = repair_report_figures(report, Path("custom.html"))
        self.assertEqual(Path(summary["html"]).resolve(), (self.outdir / "custom.html").resolve())
        self.assert_downloadable(Path(summary["html"]))

    def test_missing_image_fails_without_modifying_original_or_publishing_html(self):
        report = self.original("<img src='../05_gene_expression/missing.png'>")
        before = report.read_bytes()
        with self.assertRaisesRegex(RuntimeError, "Missing figure"):
            repair_report_figures(report)
        self.assertEqual(report.read_bytes(), before)
        self.assertFalse((self.outdir / "report.portable.html").exists())

    def test_external_path_is_rejected_and_original_cannot_be_overwritten(self):
        outside = self.root / "private.png"
        outside.write_bytes(PNG)
        report = self.original('<img src="' + outside.as_uri() + '">')
        with self.assertRaisesRegex(RuntimeError, "outside this results"):
            repair_report_figures(report)
        with self.assertRaisesRegex(RuntimeError, "separate HTML"):
            repair_report_figures(report, report)

    def test_embedded_and_remote_images_are_retained_without_fetching(self):
        body = '<img src="data:image/png;base64,AA=="><img src="https://example.test/a.png">'
        report = self.original(body)
        summary = repair_report_figures(report)
        self.assertIn(body, Path(summary["html"]).read_text())
        self.assertEqual(summary["embedded_images"], 1)
        self.assertEqual(summary["remote_images"], 1)

    def test_regenerating_report_includes_assets_in_receipt_and_repairs_missing_asset(self):
        project = {
            "project_id": "portable_test", "rseqc": {"enabled": True},
            "modules": {"gene_expression": False, "apa_a": False, "apa_a2": False,
                        "tracks": False, "rseqc": True},
        }
        reference = {"assembly": "GRCm39"}
        plan = RunPlan(project, [{"sample_id": "S1"}], [], [], reference, {"GRCm39": reference})
        with patch("rnaends2tracks.report.generate_provenance_dashboard", return_value=[]):
            report = make_report(plan, self.results)
            manifest = json.loads((self.outdir / "report_image_assets.json").read_text())
            self.assertEqual(len(manifest["images"]), 2)
            asset = self.outdir / manifest["images"][0]["asset"]
            receipt = json.loads((self.outdir / "run_receipt.json").read_text())
            self.assertIn(str(asset.resolve()), {row["path"] for row in receipt["outputs"]})
            asset.unlink()
            make_report(plan, self.results)
            self.assertTrue(asset.is_file())
        self.assert_downloadable(report)


if __name__ == "__main__":
    unittest.main()
