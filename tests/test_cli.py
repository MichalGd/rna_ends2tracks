import argparse
import csv
import json
import tempfile
import unittest
from pathlib import Path

from rnaends2tracks.cli import _downstream_branch_sequences, _status_observations, execute

HEADER = (
    "sample_id,description,genome,biological_replicate_id,technical_replicate_id,lane_id,"
    "fastq_r1,fastq_r2,condition,batch,subject,library_protocol,library_layout,read_length,"
    "kit_catalog,umi_present\n"
)


class CliTests(unittest.TestCase):
    def test_downstream_scheduler_pipelines_tracks_after_dge_and_prioritizes_apa_b(self):
        branches = _downstream_branch_sequences([
            "gene_expression", "apa_a", "apa_b", "tracks",
        ])
        self.assertEqual(branches, [
            ("apa_b", ("apa_b",)),
            ("gene_expression_then_tracks", ("gene_expression", "tracks")),
            ("apa_a", ("apa_a",)),
        ])

    def test_status_counts_enrichment_jobs_and_completed_receipts(self):
        with tempfile.TemporaryDirectory() as temporary:
            results = Path(temporary) / "results"
            index = results / "10_reports" / "enrichment_summary" / "enrichment_index.tsv"
            index.parent.mkdir(parents=True)
            with index.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=["analysis_type"], delimiter="\t")
                writer.writeheader(); writer.writerows([{"analysis_type": "dge"}, {"analysis_type": "apa_a"}])
            receipt = results / "02_alignment" / "run_receipt.json"
            receipt.parent.mkdir(parents=True)
            receipt.write_text(json.dumps({"schema_version": 1, "exit_status": 0}), encoding="utf-8")
            observations = _status_observations(results, {"workflow_status": "completed", "pid": 1})
            self.assertEqual(observations["outputs"]["enrichment"], 2)
            self.assertEqual(observations["stage_receipts"]["alignment"], "completed")

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
