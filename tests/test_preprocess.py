import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from rnaends2tracks.config import RunPlan
from rnaends2tracks.preprocess import _remove_owned_temporary_tree, preprocess


class PreprocessOrderTests(unittest.TestCase):
    def test_multiqc_read_only_temporary_tree_is_removed(self):
        with tempfile.TemporaryDirectory() as temporary:
            tree = Path(temporary) / "multiqc"
            copied_template = tree / "template"
            copied_template.mkdir(parents=True)
            (copied_template / "footer.html").write_text("footer", encoding="utf-8")
            copied_template.chmod(0o555)

            _remove_owned_temporary_tree(tree)

            self.assertFalse(tree.exists())

    def test_portable_dry_run_does_not_stat_missing_fastqs(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            results = root / "results"
            resources = {"temporary_directory": "", "preprocess": {
                "trim_parallel_jobs": 1, "star_parallel_jobs": 1, "fastqc_threads": 1,
                "bbduk_threads": 1, "bbduk_memory_gb": 1, "star_threads": 1,
                "star_memory_gb": 1, "samtools_threads": 1,
                "samtools_sort_memory_per_thread_gb": 1, "merge_parallel_jobs": 1,
                "merge_memory_gb": 1}}
            reference = {"assembly": "GRCm39", "star_index": str(root / "missing_star")}
            plan = RunPlan(
                project={"resources": resources, "preprocessing": {"minimum_length": 20, "trim_quality": 10},
                         "protocol": {"orientation_min_fraction": 0.75}},
                samples=[{"sample_id": "S1", "genome": "GRCm39"}],
                sample_rows=[{"sample_id": "S1", "genome": "GRCm39", "technical_replicate_id": "T01",
                              "lane_id": "L001", "fastq_r1": str(root / "missing.fastq.gz")}],
                contrasts=[], reference=reference, references={"GRCm39": reference})
            preprocess(plan, results, dry_run=True)
            self.assertIn("DRY RUN: STAR", (results / "logs" / "preprocess" / "S1.T01.L001.star.log").read_text())

    def test_trim_and_star_use_separate_bounded_phases(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            fastq = root / "sample.fastq.gz"; fastq.write_bytes(b"test-fastq")
            results = root / "results"
            resources = {"temporary_directory": "", "preprocess": {
                "trim_parallel_jobs": 2, "star_parallel_jobs": 1, "fastqc_threads": 1,
                "bbduk_threads": 1, "bbduk_memory_gb": 1, "star_threads": 1,
                "star_memory_gb": 1, "samtools_threads": 1,
                "samtools_sort_memory_per_thread_gb": 1, "merge_parallel_jobs": 1,
                "merge_memory_gb": 1}}
            reference = {
                "assembly": "GRCm39", "star_index": str(root / "star"),
                "fasta": str(root / "genome.fa"), "gtf": str(root / "genes.gtf"),
                "chrom_sizes": str(root / "chrom.sizes"),
            }
            plan = RunPlan(
                project={"resources": resources, "preprocessing": {"minimum_length": 20, "trim_quality": 10},
                         "protocol": {"orientation_min_fraction": 0.75}},
                samples=[{"sample_id": "S1", "genome": "GRCm39"}],
                sample_rows=[{"sample_id": "S1", "genome": "GRCm39", "technical_replicate_id": "T01",
                              "lane_id": "L001", "fastq_r1": str(fastq)}],
                contrasts=[], reference=reference, references={"GRCm39": reference})
            calls, commands, environments, phases = [], [], [], []

            def fake_run(command, _log, dry_run=False, cwd=None, env=None):
                self.assertFalse(dry_run); calls.append(command[0]); commands.append(command); environments.append(env)
                if command[0] == "bbduk.sh":
                    Path(next(value[4:] for value in command if value.startswith("out="))).write_bytes(b"trimmed")
                elif command[0] == "STAR":
                    prefix = command[command.index("--outFileNamePrefix") + 1]
                    Path(prefix + "ReadsPerGene.out.tab").write_text("gene1\t100\t10\t90\n", encoding="utf-8")
                    Path(prefix + "Aligned.out.bam").write_bytes(b"unsorted")
                elif command[:2] in (["samtools", "sort"], ["samtools", "view"]):
                    Path(command[command.index("-o") + 1]).write_bytes(b"bam")
                elif command[:2] == ["samtools", "index"]:
                    Path(command[command.index("-o") + 1]).write_bytes(b"index")

            def immediate(stage, jobs, workers, _timing):
                phases.append((stage, workers)); return [worker() for _, worker in jobs]

            def fake_run_to_path(_command, output, _log, dry_run=False, cwd=None, env=None):
                self.assertFalse(dry_run)
                Path(output).write_text("flagstat\n", encoding="utf-8")

            with (patch("rnaends2tracks.preprocess.require_tools"),
                  patch("rnaends2tracks.preprocess.run", side_effect=fake_run),
                  patch("rnaends2tracks.preprocess.run_to_path", side_effect=fake_run_to_path),
                  patch("rnaends2tracks.preprocess.run_bounded", side_effect=immediate),
                  patch("rnaends2tracks.preprocess.write_receipt"),
                  patch("rnaends2tracks.preprocess.signature_for", return_value="signature")):
                preprocess(plan, results)

            self.assertEqual(phases[:2], [("qc_and_trim", 2), ("star_and_sort", 1)])
            self.assertLess(calls.index("bbduk.sh"), calls.index("STAR"))
            bbduk = next(command for command in commands if command[0] == "bbduk.sh")
            self.assertIn("qtrim=r", bbduk); self.assertNotIn("qtrim=t", bbduk)
            orientation = (results / "02_alignment" / "protocol_orientation.tsv").read_text(encoding="utf-8")
            self.assertIn("S1\tT01\tL001\t10\t90\t0.9\tpass", orientation)
            star = next(command for command in commands if command[0] == "STAR")
            self.assertIn("ID:S1.T01.L001", star); self.assertIn("LB:S1.T01", star)
            multiqc_index = calls.index("multiqc")
            self.assertIn("--no-clean-up", commands[multiqc_index])
            self.assertIn(str(results / ".checkpoints" / "multiqc_tmp"),
                          environments[multiqc_index]["TMPDIR"])
            self.assertFalse((results / ".checkpoints" / "multiqc_tmp").joinpath(
                Path(environments[multiqc_index]["TMPDIR"]).name).exists())
            self.assertTrue((results / "02_alignment" / "S1" / "S1.bam").is_file())
            self.assertTrue((results / "02_alignment" / "S1" / "S1.bam.bai").is_file())


if __name__ == "__main__":
    unittest.main()
