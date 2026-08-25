import csv
import gzip
import tempfile
import unittest
from pathlib import Path

import yaml

from rnaends2tracks.config import (
    ConfigError,
    REQUIRED_COLUMNS,
    build_plan,
    collapse_samples,
    generate_contrasts,
    load_samplesheet,
    resolve_contrast_designs,
    validate_design,
)


def sample(sample_id, condition, batch, subject=""):
    return {
        "sample_id": sample_id, "description": f"Description for {sample_id}", "genome": "GRCh38",
        "biological_replicate_id": sample_id, "technical_replicate_id": "T01", "condition": condition,
        "batch": batch, "subject": subject, "library_protocol": "quantseq_rev_v2_se",
        "library_layout": "SE", "read_length": "100", "kit_catalog": "REV_V2", "umi_present": "false",
    }


class ContractTests(unittest.TestCase):
    def test_repository_example_resolves_complete_pairs(self):
        root = Path(__file__).resolve().parents[1]
        plan = build_plan(
            root / "config" / "project.example.yaml",
            root / "config" / "samplesheet.example.csv",
            check_inputs=False,
        )
        self.assertEqual(len(plan.contrasts), 1)
        self.assertEqual(plan.contrasts[0]["design_mode"], "paired")
        self.assertEqual(plan.contrasts[0]["n_pairs"], 2)
        self.assertEqual(plan.contrasts[0]["resolved_design"], "~ subject + condition")

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

    def test_one_level_design_term_fails_before_r(self):
        samples = [sample("A1", "A", "run8"), sample("A2", "A", "run8"),
                   sample("B1", "B", "run8"), sample("B2", "B", "run8")]
        with self.assertRaisesRegex(ConfigError, "fewer than two levels: batch"):
            validate_design(samples, "~ batch + condition")

    def test_nested_subjects_resolve_three_paired_and_twelve_unpaired_contrasts(self):
        samples = []
        order = []
        for target in ("CTCF", "RAD21", "WAPL"):
            control = f"{target}_control"; treated = f"{target}_IAA"
            order.extend((control, treated))
            for replicate in (1, 2, 3):
                subject = f"{target}_R{replicate}"
                samples.append(sample(f"{control}_R{replicate}", control, "run8", subject))
                samples.append(sample(f"{treated}_R{replicate}", treated, "run8", subject))
        project = {
            "design": "~ condition",
            "statistics": {"pairing": {
                "mode": "auto", "subject_column": "subject",
                "paired_design": "~ subject + condition", "incomplete_pair_action": "error",
            }},
        }
        contrasts = resolve_contrast_designs(samples, generate_contrasts(samples, order), project)
        paired = [contrast for contrast in contrasts if contrast["design_mode"] == "paired"]
        unpaired = [contrast for contrast in contrasts if contrast["design_mode"] == "unpaired"]
        self.assertEqual(len(contrasts), 15)
        self.assertEqual(
            [contrast["contrast_id"] for contrast in paired],
            ["CTCF_IAA_vs_CTCF_control", "RAD21_IAA_vs_RAD21_control", "WAPL_IAA_vs_WAPL_control"],
        )
        self.assertTrue(all(contrast["resolved_design"] == "~ subject + condition" for contrast in paired))
        self.assertTrue(all(contrast["n_pairs"] == 3 for contrast in paired))
        self.assertEqual(len(unpaired), 12)
        self.assertTrue(all(contrast["resolved_design"] == "~ condition" for contrast in unpaired))
        self.assertTrue(all(contrast["pairing_status"] == "disjoint_subjects" for contrast in unpaired))

    def test_incomplete_pairing_fails_in_auto_error_mode(self):
        samples = [
            sample("A1", "A", "run8", "S1"), sample("A2", "A", "run8", "S2"),
            sample("B1", "B", "run8", "S1"), sample("B2", "B", "run8", "S3"),
        ]
        project = {
            "design": "~ condition",
            "statistics": {"pairing": {"mode": "auto", "incomplete_pair_action": "error"}},
        }
        with self.assertRaisesRegex(ConfigError, "Incomplete subject pairing"):
            resolve_contrast_designs(samples, generate_contrasts(samples, ["A", "B"]), project)

    def test_pairing_none_preserves_legacy_unpaired_design(self):
        samples = [
            sample("A1", "A", "run8", "S1"), sample("A2", "A", "run8", "S2"),
            sample("B1", "B", "run8", "S1"), sample("B2", "B", "run8", "S2"),
        ]
        project = {"design": "~ condition", "statistics": {"pairing": {"mode": "none"}}}
        contrast = resolve_contrast_designs(samples, generate_contrasts(samples, ["A", "B"]), project)[0]
        self.assertTrue(contrast["paired"])
        self.assertEqual(contrast["design_mode"], "unpaired")
        self.assertEqual(contrast["resolved_design"], "~ condition")

    def test_pairing_column_must_be_a_valid_formula_term(self):
        samples = [sample("A1", "A", "X"), sample("A2", "A", "Y"),
                   sample("B1", "B", "X"), sample("B2", "B", "Y")]
        project = {"design": "~ condition", "statistics": {"pairing": {
            "mode": "auto", "subject_column": "subject-id",
        }}}
        with self.assertRaisesRegex(ConfigError, "safe column name"):
            resolve_contrast_designs(samples, generate_contrasts(samples, ["A", "B"]), project)

    def test_technical_replicates_and_lanes_collapse_to_one_biological_sample(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            with (root / "samples.csv").open("w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(handle, fieldnames=list(REQUIRED_COLUMNS)); writer.writeheader()
                for technical, lane in (("T01", "L001"), ("T01", "L002"), ("T02", "L001")):
                    row = sample("S1", "A", "B1")
                    row.update({
                        "technical_replicate_id": technical, "lane_id": lane,
                        "fastq_r1": f"S1_{technical}_{lane}.fastq.gz", "fastq_r2": "",
                    })
                    writer.writerow(row)
            rows = load_samplesheet(root / "samples.csv", check_fastqs=False)
            collapsed = collapse_samples(rows)
            self.assertEqual(len(collapsed), 1)
            self.assertEqual(collapsed[0]["technical_replicate_count"], "2")
            self.assertEqual(collapsed[0]["sequencing_lane_count"], "3")
            self.assertEqual(collapsed[0]["description"], "Description for S1")
            self.assertNotIn("technical_replicate_id", collapsed[0])

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
            fields = list(REQUIRED_COLUMNS)
            with (root / "samples.csv").open("w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(handle, fieldnames=fields); writer.writeheader()
                for condition in ("A", "B"):
                    for index in (1, 2):
                        row = sample(f"{condition}{index}", condition, "B1")
                        row.update({
                            "genome": "mm39", "lane_id": "L001",
                            "fastq_r1": f"{condition}{index}.fastq.gz", "fastq_r2": "",
                        })
                        writer.writerow(row)
            plan = build_plan(root / "project.yaml", root / "samples.csv", check_inputs=False)
            self.assertEqual(plan.reference["species"], "mouse")
            self.assertEqual({item["genome"] for item in plan.samples}, {"GRCm39"})
            self.assertEqual(len(plan.contrasts), 1)
            original = (root / "samples.csv").read_text(encoding="utf-8")
            mixed = original.replace("mm39", "GRCh38", 1)
            (root / "samples.csv").write_text(mixed, encoding="utf-8")
            with self.assertRaisesRegex(ConfigError, "exactly one genome"):
                build_plan(root / "project.yaml", root / "samples.csv", check_inputs=False)
            (root / "samples.csv").write_text(original.replace("mm39", "GRCh38"), encoding="utf-8")
            with self.assertRaisesRegex(ConfigError, "does not match reference-manifest assembly"):
                build_plan(root / "project.yaml", root / "samples.csv", check_inputs=False)

    def test_existing_star_index_is_reused_only_when_contigs_and_lengths_match(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary); star = root / "star"; star.mkdir()
            (root / "genome.fa").write_text(">chr1\n" + "A" * 100 + "\n", encoding="utf-8")
            (root / "genome.fa.fai").write_text("chr1\t100\t6\t100\t101\n", encoding="utf-8")
            (root / "genes.gtf").write_text('chr1\ttest\tgene\t1\t100\t.\t+\t.\tgene_id "g1";\n', encoding="utf-8")
            (root / "sizes").write_text("chr1\t100\n", encoding="utf-8")
            for asset in ("Genome", "SA", "SAindex"):
                (star / asset).write_bytes(b"index")
            (star / "chrName.txt").write_text("chr1\n", encoding="utf-8")
            (star / "chrLength.txt").write_text("100\n", encoding="utf-8")
            (star / "genomeParameters.txt").write_text("versionGenome\t2.7.11b\nsjdbOverhang\t99\n", encoding="utf-8")
            reference = {"species": "human", "assembly": "GRCh38", "fasta": "genome.fa", "gtf": "genes.gtf",
                         "star_index": "star", "chrom_sizes": "sizes"}
            (root / "reference.yaml").write_text(yaml.safe_dump(reference), encoding="utf-8")
            project = {"project_id": "star_reuse", "condition_order": ["A", "B"], "design": "~ condition",
                       "reference": {"manifest": "reference.yaml"},
                       "protocol": {"profile": "quantseq_rev_v2_se", "has_umi": False, "retain_duplicate_flagged_reads": True}}
            (root / "project.yaml").write_text(yaml.safe_dump(project), encoding="utf-8")
            fields = list(REQUIRED_COLUMNS)
            with (root / "samples.csv").open("w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(handle, fieldnames=fields); writer.writeheader()
                for condition in ("A", "B"):
                    for index in (1, 2):
                        fastq = root / f"{condition}{index}.fastq.gz"
                        with gzip.open(fastq, "wt", encoding="ascii") as fq:
                            fq.write("@r1\n" + "A" * 20 + "\n+\n" + "I" * 20 + "\n")
                        row = sample(f"{condition}{index}", condition, "B1")
                        row.update({"lane_id": "L001", "fastq_r1": str(fastq), "fastq_r2": ""})
                        writer.writerow(row)
            plan = build_plan(root / "project.yaml", root / "samples.csv", check_inputs=True)
            self.assertEqual(plan.reference["star_index_version"], "2.7.11b")
            self.assertEqual(plan.reference["star_sjdb_overhang"], 99)
            self.assertEqual(plan.reference["_warnings"][0]["warning_code"], "STAR_INDEX_PROVENANCE_UNVERIFIED")
            (star / "chrLength.txt").write_text("101\n", encoding="utf-8")
            with self.assertRaises(ConfigError):
                build_plan(root / "project.yaml", root / "samples.csv", check_inputs=True)


if __name__ == "__main__":
    unittest.main()
