from __future__ import annotations

import csv
import gzip
import json
import os
import subprocess
import sys
import tempfile
import unittest
import zipfile
from pathlib import Path

from rnaends2tracks.pas_sources import (
    PasSourceRecord,
    lift_records,
    prepare_sources,
    read_gencode_polya_sites,
    read_polyadb_collection,
)


class PasSourceTests(unittest.TestCase):
    @staticmethod
    def _write_sources(root: Path) -> tuple[Path, Path]:
        gencode = root / "polyAs.gtf.gz"
        with gzip.open(gencode, "wt", encoding="utf-8") as handle:
            handle.write('chr1\tHAVANA\tpolyA_site\t100\t101\t.\t+\t.\tgene_id "plus";\n')
        polyadb = root / "HumanPas.zip"
        with zipfile.ZipFile(polyadb, "w") as archive:
            archive.writestr(
                "HumanPas/hg38.PAS.main.tsv",
                "PAS_ID\tGeneSymbol\nchr1:+:101\tGeneA\n",
            )
            archive.writestr(
                "HumanPas/hg38.PAS.max.tsv",
                "PAS_ID\tAvgRPM\tPSE\tGeneSymbol\nchr1:+:201\t1\t2\tGeneB\n",
            )
        return gencode, polyadb

    def test_gencode_retains_only_cleavage_sites_with_strand_aware_coordinates(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "polyAs.gtf.gz"
            with gzip.open(path, "wt", encoding="utf-8") as handle:
                handle.write('chr1\tHAVANA\tpolyA_signal\t90\t95\t.\t+\t.\tgene_id "signal";\n')
                handle.write('chr1\tHAVANA\tpolyA_site\t100\t101\t.\t+\t.\tgene_id "plus";\n')
                handle.write('chr2\tHAVANA\tpolyA_site\t200\t201\t.\t-\t.\tgene_id "minus";\n')
                handle.write('chr2\tHAVANA\tpseudo_polyA\t300\t301\t.\t-\t.\tgene_id "pseudo";\n')
            records, audit = read_gencode_polya_sites(path)
            self.assertEqual([(r.chrom, r.position, r.strand) for r in records],
                             [("chr1", 100, "+"), ("chr2", 199, "-")])
            self.assertEqual(audit["retained_unique_polyA_site"], 2)
            self.assertEqual(audit["feature_polyA_signal"], 1)
            self.assertEqual(audit["feature_pseudo_polyA"], 1)

    def test_polyadb_collections_parse_one_based_ids_and_deduplicate_gene_rows(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "MousePas.zip"
            with zipfile.ZipFile(path, "w") as archive:
                archive.writestr(
                    "MousePas/mm10.PAS.main.tsv",
                    "PAS_ID\tGeneSymbol\nchr1:-:11\tGeneA\nchr1:-:11\tGeneB\nchr2:+:21\tGeneC\n",
                )
                archive.writestr(
                    "MousePas/mm10.PAS.max.tsv",
                    "PAS_ID\tAvgRPM\tPSE\tGeneSymbol\nchr3:+:31\t1\t2\tGeneD\n",
                )
            main, audit = read_polyadb_collection(path, "main")
            maximum, _ = read_polyadb_collection(path, "max")
            self.assertEqual([(r.chrom, r.position, r.strand) for r in main],
                             [("chr1", 10, "-"), ("chr2", 20, "+")])
            self.assertEqual([(r.chrom, r.position, r.strand) for r in maximum],
                             [("chr3", 30, "+")])
            self.assertEqual(audit["duplicate_gene_assignment_rows"], 1)

    def test_prepare_human_sources_writes_builder_manifest_and_deterministic_bed(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            gencode, polyadb = self._write_sources(root)
            output = root / "prepared"
            provenance = prepare_sources(
                species="human",
                assembly="GRCh38",
                annotation_release="GENCODE_v42",
                gencode_polya_gtf=gencode,
                polyadb_zip=polyadb,
                polyadb_assembly="hg38",
                output=output,
                download_date="2026-08-26",
            )
            manifest = json.loads((output / "source_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual(len(manifest["sources"]), 3)
            self.assertEqual(provenance["gencode_feature_policy"], "polyA_site_only")
            with gzip.open(output / "gencode_polyA.GRCh38.bed.gz", "rt", encoding="utf-8") as handle:
                self.assertEqual(handle.readline().split("\t")[:6],
                                 ["chr1", "100", "101", "plus", "0", "+"])

    def test_mouse_polyadb_requires_target_liftover_chain(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            gencode, polyadb = self._write_sources(root)
            with self.assertRaisesRegex(ValueError, "requires a lift-over chain"):
                prepare_sources(
                    species="mouse",
                    assembly="GRCm39",
                    annotation_release="GENCODE_vM31",
                    gencode_polya_gtf=gencode,
                    polyadb_zip=polyadb,
                    polyadb_assembly="mm10",
                    output=root / "prepared_mouse",
                    download_date="2026-08-26",
                )

    @unittest.skipIf(os.name == "nt", "POSIX executable fixture")
    def test_liftover_keeps_unique_strand_preserved_sites(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            chain = root / "test.chain"
            chain.write_text("fixture\n", encoding="utf-8")
            executable = root / "liftOver"
            executable.write_text(
                "#!/usr/bin/env bash\nset -euo pipefail\ncp \"$2\" \"$4\"\n: > \"$5\"\n",
                encoding="utf-8",
            )
            executable.chmod(0o755)
            records, audit = lift_records(
                [PasSourceRecord("chr1", 10, "+", "PAS1")], chain, str(executable)
            )
            self.assertEqual([(row.chrom, row.position, row.strand) for row in records],
                             [("chr1", 10, "+")])
            self.assertEqual(audit["unique_strand_preserved"], 1)

    def test_prepared_sources_build_core_and_rescue_atlas(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            gencode, polyadb = self._write_sources(root)
            prepared = root / "prepared"
            prepare_sources(
                species="human",
                assembly="GRCh38",
                annotation_release="GENCODE_v42",
                gencode_polya_gtf=gencode,
                polyadb_zip=polyadb,
                polyadb_assembly="hg38",
                output=prepared,
                download_date="2026-08-26",
            )
            target_gtf = root / "target.gtf"
            target_gtf.write_text(
                'chr1\tGENCODE\tgene\t1\t1000\t.\t+\t.\tgene_id "GENE1";\n'
                'chr1\tGENCODE\texon\t1\t1000\t.\t+\t.\tgene_id "GENE1";\n',
                encoding="utf-8",
            )
            chrom_sizes = root / "chrom.sizes"
            chrom_sizes.write_text("chr1\t1000\n", encoding="utf-8")
            atlas = root / "atlas"
            project = Path(__file__).resolve().parents[1]
            bootstrap = (
                "import runpy,sys,types; "
                "module=types.ModuleType('pysam'); module.AlignedSegment=object; "
                "sys.modules['pysam']=module; yaml=types.ModuleType('yaml'); "
                "yaml.safe_load=lambda value: {}; sys.modules['yaml']=yaml; "
                "script=sys.argv[1]; sys.argv=sys.argv[1:]; "
                "runpy.run_path(script,run_name='__main__')"
            )
            subprocess.run(
                [
                    sys.executable, "-c", bootstrap,
                    str(project / "scripts" / "build_pas_atlas.py"),
                    "--species", "human", "--assembly", "GRCh38",
                    "--annotation-release", "GENCODE_v42",
                    "--atlas-id", "test_GRCh38_v1",
                    "--source-manifest", str(prepared / "source_manifest.json"),
                    "--gtf", str(target_gtf), "--chrom-sizes", str(chrom_sizes),
                    "--gencode", str(prepared / "gencode_polyA.GRCh38.bed.gz"),
                    "--polyadb-main", str(prepared / "polyadb_v4.1_main.GRCh38.bed.gz"),
                    "--polyadb-max", str(prepared / "polyadb_v4.1_max.GRCh38.bed.gz"),
                    "--output", str(atlas),
                ],
                check=True,
            )
            with gzip.open(atlas / "master.tsv.gz", "rt", encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual([(row["tier"], row["confidence"]) for row in rows],
                             [("core", "A"), ("rescue", "C")])
            self.assertTrue((atlas / "SHA256SUMS").is_file())


if __name__ == "__main__":
    unittest.main()
