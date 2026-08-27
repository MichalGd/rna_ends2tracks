import argparse
import tempfile
import unittest
from pathlib import Path

from rnaends2tracks.cli import execute

HEADER = (
    "sample_id,description,genome,biological_replicate_id,technical_replicate_id,lane_id,"
    "fastq_r1,fastq_r2,condition,batch,subject,library_protocol,library_layout,read_length,"
    "kit_catalog,umi_present\n"
)


class CliTests(unittest.TestCase):
    def test_complete_portable_dry_run_accepts_missing_inputs(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            rows = []
            for condition in ("control", "treated"):
                for replicate in (1, 2):
                    sample = f"{condition}_{replicate}"
                    rows.append(
                        f"{sample},example,GRCm39,{sample},T01,L001,{root / (sample + '.fastq.gz')},,"
                        f"{condition},B1,{sample},quantseq_rev_v2_se,SE,101,REV_V2,false"
                    )
            sheet = root / "samplesheet.csv"
            sheet.write_text(HEADER + "\n".join(rows) + "\n", encoding="utf-8")
            config = root / "config.conf"
            config.write_text("\n".join([
                "PROJECT_ID=portable_dry_run",
                f'SAMPLESHEET="{sheet}"',
                f'OUTPUT_DIR="{root / "results"}"',
                "MM39_STAR_INDEX=/missing/star",
                "MM39_FASTA=/missing/genome.fa",
                "MM39_GTF=/missing/genes.gtf",
                "MM39_CHROM_SIZES=/missing/chrom.sizes",
                "MM39_PAS_ATLAS=/missing/pas_atlas",
            ]) + "\n", encoding="utf-8")
            status = execute(argparse.Namespace(
                config=str(config), dry_run=True, from_step=None, stop_after=None,
                force_step=[], skip_input_checks=True,
            ))
            self.assertEqual(status, 0)
            self.assertTrue((root / "results" / "00_metadata" / "resource_plan.tsv").is_file())
            events = (root / "results" / "logs" / "events.jsonl").read_text(encoding="utf-8")
            self.assertIn('"module": "preprocess"', events)
            self.assertIn('"module": "c0_tracks"', events)
            self.assertIn('"status": "dry_run"', events)
            master = (root / "results" / "rna_ends2tracks.log").read_text(encoding="utf-8")
            self.assertIn("[workflow] STARTED", master)
            self.assertIn("[workflow] COMPLETED Dry run completed", master)
            status = (root / "results" / "00_metadata" / "run_status.json").read_text(encoding="utf-8")
            self.assertIn('"workflow_status": "completed"', status)


if __name__ == "__main__":
    unittest.main()
