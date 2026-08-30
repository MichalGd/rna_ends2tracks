import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from rnaends2tracks.config import RunPlan
from rnaends2tracks.preprocess import (
    _c0_overlap_workers, _remove_owned_temporary_tree, _run_fastq_screen, preprocess,
)


class PreprocessOrderTests(unittest.TestCase):
    def test_fastq_screen_processes_both_paired_mates_and_records_reports(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary); results = root / "results"
            config = root / "fastq_screen.conf"; config.write_text("DATABASE\tMouse\t/index/mm39\n")
            row = {
                "sample_id": "S1", "genome": "GRCm39", "technical_replicate_id": "T01",
                "lane_id": "L001", "fastq_r1": str(root / "S1_R1.fastq.gz"),
                "fastq_r2": str(root / "S1_R2.fastq.gz"), "library_layout": "PE",
            }
            for key in ("fastq_r1", "fastq_r2"):
                Path(row[key]).write_bytes(b"fastq")
            plan = RunPlan(
                project={
                    "preprocessing": {"fastq_screen": {
                        "enabled": True, "config": str(config), "missing_action": "error", "subset": 1000,
                    }},
                    "resources": {"preprocess": {
                        "fastq_screen_parallel_jobs": 1, "fastq_screen_threads": 2,
                        "fastq_screen_memory_gb": 4,
                    }},
                },
                samples=[{"sample_id": "S1", "genome": "GRCm39", "library_layout": "PE"}],
                sample_rows=[row], contrasts=[], reference={"assembly": "GRCm39"},
                references={"GRCm39": {"assembly": "GRCm39"}},
            )
            commands = []

            def fake_run(command, _log, dry_run=False, cwd=None, env=None):
                commands.append(command)
                outdir = Path(command[command.index("--outdir") + 1])
                path = command[-1]
                (outdir / f"{Path(path).name.removesuffix('.gz').removesuffix('.fastq')}_screen.txt").write_text(
                    "FastQ Screen report\nGenome\t#Reads_processed\t%Unmapped\t%One_hit_one_genome\t"
                    "%Multiple_hits_one_genome\t%One_hit_multiple_genomes\t"
                    "%Multiple_hits_multiple_genomes\nMouse\t1000\t5\t90\t2\t2\t1\n",
                    encoding="utf-8",
                )

            def immediate(_stage, jobs, _workers, _timing, progress=None):
                return [worker() for _label, worker in jobs]

            with (patch("rnaends2tracks.preprocess.require_tools"),
                  patch("rnaends2tracks.preprocess.run", side_effect=fake_run),
                  patch("rnaends2tracks.preprocess.run_bounded", side_effect=immediate),
                  patch("rnaends2tracks.preprocess.signature_for", return_value="signature"),
                  patch("rnaends2tracks.preprocess.write_receipt")):
                summary = _run_fastq_screen(plan, results, False, False, None)

            self.assertEqual(len(commands), 2)
            self.assertEqual([command[-1] for command in commands], [row["fastq_r1"], row["fastq_r2"]])
            self.assertEqual([Path(command[command.index("--outdir") + 1]).name for command in commands],
                             ["R1", "R2"])
            text = summary.read_text(encoding="utf-8")
            self.assertIn("\tPE\tR1,R2\tPASS\t", text)
            self.assertEqual(text.count("_screen.txt"), 2)
            metrics = summary.with_name("fastq_screen_metrics.tsv").read_text(encoding="utf-8")
            self.assertIn("\tR1\tMouse\t1000\t5\t90\t2\t2\t1\t", metrics)
            self.assertIn("\tR2\tMouse\t1000\t5\t90\t2\t2\t1\t", metrics)

    def test_paired_dry_run_uses_mate_aware_bbduk_and_star(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary); results = root / "results"
            resources = {"temporary_directory": "", "preprocess": {
                "trim_parallel_jobs": 1, "star_parallel_jobs": 1, "fastqc_threads": 1,
                "bbduk_threads": 1, "bbduk_memory_gb": 1, "star_threads": 1,
                "star_memory_gb": 1, "samtools_threads": 1,
                "samtools_sort_memory_per_thread_gb": 1, "merge_parallel_jobs": 1,
                "merge_memory_gb": 1}}
            reference = {"assembly": "GRCm39", "star_index": str(root / "star")}
            row = {"sample_id": "S1", "genome": "GRCm39", "technical_replicate_id": "T01",
                   "lane_id": "L001", "fastq_r1": str(root / "R1.fastq.gz"),
                   "fastq_r2": str(root / "R2.fastq.gz"), "library_layout": "PE"}
            plan = RunPlan(
                project={"resources": resources, "preprocessing": {"minimum_length": 20, "trim_quality": 10},
                         "protocol": {"orientation_min_fraction": 0.75, "library_layout": "PE"}},
                samples=[{"sample_id": "S1", "genome": "GRCm39", "library_layout": "PE"}],
                sample_rows=[row], contrasts=[], reference=reference, references={"GRCm39": reference})
            preprocess(plan, results, dry_run=True)
            trim_log = (results / "logs" / "preprocess" / "S1.T01.L001.trim.log").read_text()
            star_log = (results / "logs" / "preprocess" / "S1.T01.L001.star.log").read_text()
            self.assertIn("in1=", trim_log); self.assertIn("in2=", trim_log)
            self.assertIn("out1=", trim_log); self.assertIn("out2=", trim_log)
            self.assertIn("ftl=12", trim_log); self.assertIn("skipr1=t", trim_log)
            self.assertIn("R1.trimmed.fastq.gz", star_log)
            self.assertIn("R2.trimmed.fastq.gz", star_log)

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

            def immediate(stage, jobs, workers, _timing, progress=None, on_completed=None):
                phases.append((stage, workers))
                values = []
                for label, worker in jobs:
                    value = worker()
                    values.append(value)
                    if on_completed is not None:
                        on_completed(label, value)
                    if progress is not None:
                        progress(label, "completed")
                return values

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

    def test_c0_overlap_workers_respect_combined_cpu_and_memory_budgets(self):
        samples = [
            {"sample_id": "S1", "genome": "GRCm39"},
            {"sample_id": "S2", "genome": "GRCm39"},
        ]
        resources = {
            "total_threads": 16,
            "total_memory_gb": 32,
            "preprocess": {
                "merge_parallel_jobs": 2, "samtools_threads": 2, "merge_memory_gb": 4,
            },
            "tracks": {"parallel_jobs": 4, "samtools_threads": 2, "memory_gb": 4},
        }
        project = {
            "modules": {"tracks": True},
            "resources": resources,
            "tracks": {"early_c0": True, "families": {"all_reads": True}},
        }
        reference = {"assembly": "GRCm39"}
        plan = RunPlan(project, samples, [], [], reference, {"GRCm39": reference})

        self.assertEqual(_c0_overlap_workers(plan), 2)
        resources["total_threads"] = 4
        self.assertEqual(_c0_overlap_workers(plan), 0)


if __name__ == "__main__":
    unittest.main()
