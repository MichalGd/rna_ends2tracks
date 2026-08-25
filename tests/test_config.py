import csv
import tempfile
import unittest
from pathlib import Path

import yaml

from rnaends2tracks.config import ConfigError, build_plan, generate_contrasts, validate_design


def sample(sample_id, condition, batch, subject=""):
    return {
        "sample_id": sample_id, "biological_replicate_id": sample_id, "condition": condition,
        "batch": batch, "subject": subject, "library_protocol": "quantseq_rev_v2_se",
        "library_layout": "SE", "read_length": "100", "kit_catalog": "REV_V2", "umi_present": "false",
    }


class ContractTests(unittest.TestCase):
    def test_all_pairs_and_n2_warning(self):
        samples = [sample(f"{condition}_{i}", condition, f"B{i}") for condition in ("A", "B", "C") for i in (1, 2)]
        contrasts = generate_contrasts(samples, ["A", "B", "C"])
        self.assertEqual([x["contrast_id"] for x in contrasts], ["B_vs_A", "C_vs_A", "C_vs_B"])
        self.assertTrue(all(x["design_status"] == "LOW_REPLICATION_N2" for x in contrasts))

    def test_confounded_design_fails(self):
        samples = [sample("A1", "A", "X"), sample("A2", "A", "X"), sample("B1", "B", "Y"), sample("B2", "B", "Y")]
        with self.assertRaises(ConfigError):
            validate_design(samples, "~ batch + condition")

    def test_balanced_design_passes(self):
        samples = [sample("A1", "A", "X"), sample("A2", "A", "Y"), sample("B1", "B", "X"), sample("B2", "B", "Y")]
        validate_design(samples, "~ batch + condition")

    def test_mouse_project_and_no_umi_contract(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary); (root / "star").mkdir()
            for name in ("genome.fa", "genes.gtf", "sizes"):
                (root / name).write_text("x\n", encoding="utf-8")
            reference = {"species": "mouse", "assembly": "GRCm39", "fasta": "genome.fa", "gtf": "genes.gtf",
                         "star_index": "star", "chrom_sizes": "sizes"}
            (root / "reference.yaml").write_text(yaml.safe_dump(reference), encoding="utf-8")
            project = {"project_id": "mouse_test", "condition_order": ["A", "B"], "design": "~ condition",
                       "reference": {"manifest": "reference.yaml"},
                       "protocol": {"profile": "quantseq_rev_v2_se", "has_umi": False, "retain_duplicate_flagged_reads": True}}
            (root / "project.yaml").write_text(yaml.safe_dump(project), encoding="utf-8")
            fields = ["sample_id", "biological_replicate_id", "lane_id", "fastq_r1", "fastq_r2", "condition", "batch",
                      "subject", "library_protocol", "library_layout", "read_length", "kit_catalog", "umi_present"]
            with (root / "samples.csv").open("w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(handle, fieldnames=fields); writer.writeheader()
                for condition in ("A", "B"):
                    for index in (1, 2):
                        row = sample(f"{condition}{index}", condition, "B1")
                        row.update({"lane_id": "L001", "fastq_r1": f"{condition}{index}.fastq.gz", "fastq_r2": ""})
                        writer.writerow(row)
            plan = build_plan(root / "project.yaml", root / "samples.csv", check_inputs=False)
            self.assertEqual(plan.reference["species"], "mouse")
            self.assertEqual(len(plan.contrasts), 1)


if __name__ == "__main__":
    unittest.main()
