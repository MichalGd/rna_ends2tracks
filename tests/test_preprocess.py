import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from rnaends2tracks.config import RunPlan
from rnaends2tracks.preprocess import preprocess


class PreprocessOrderTests(unittest.TestCase):
    def test_non_dry_run_reads_orientation_counts_after_star(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            fastq = root / "sample.fastq.gz"
            fastq.write_bytes(b"test-fastq")
            results = root / "results"
            plan = RunPlan(
                project={
                    "resources": {"threads": 1},
                    "preprocessing": {"minimum_length": 20, "trim_quality": 10},
                    "protocol": {"orientation_min_fraction": 0.75},
                },
                samples=[{"sample_id": "S1"}],
                sample_rows=[{"sample_id": "S1", "lane_id": "L001", "fastq_r1": str(fastq)}],
                contrasts=[],
                reference={"star_index": str(root / "star")},
            )
            calls = []
            commands = []

            def fake_run(command, _log, dry_run=False, cwd=None):
                self.assertFalse(dry_run)
                self.assertIsNone(cwd)
                calls.append(command[0])
                commands.append(command)
                if command[0] == "bbduk.sh":
                    output = next(value[4:] for value in command if value.startswith("out="))
                    Path(output).write_bytes(b"trimmed")
                elif command[0] == "STAR":
                    prefix = command[command.index("--outFileNamePrefix") + 1]
                    Path(prefix + "ReadsPerGene.out.tab").write_text(
                        "gene1\t100\t10\t90\n", encoding="utf-8"
                    )
                    Path(prefix + "Aligned.out.bam").write_bytes(b"unsorted")
                elif command[:2] == ["samtools", "sort"]:
                    Path(command[command.index("-o") + 1]).write_bytes(b"sorted")
                elif command[:2] == ["samtools", "index"]:
                    Path(command[-1] + ".bai").write_bytes(b"index")

            with (
                patch("rnaends2tracks.preprocess.require_tools"),
                patch("rnaends2tracks.preprocess.run", side_effect=fake_run),
                patch("rnaends2tracks.preprocess.write_receipt"),
                patch("subprocess.run"),
            ):
                preprocess(plan, results)

            self.assertLess(calls.index("STAR"), calls.index("samtools"))
            bbduk_command = next(command for command in commands if command[0] == "bbduk.sh")
            self.assertIn("qtrim=r", bbduk_command)
            self.assertNotIn("qtrim=t", bbduk_command)
            orientation = (results / "02_alignment" / "protocol_orientation.tsv").read_text(encoding="utf-8")
            self.assertIn("S1\tL001\t10\t90\t0.9\tpass", orientation)
            self.assertTrue((results / "02_alignment" / "S1" / "S1.bam").is_file())
            self.assertTrue((results / "02_alignment" / "S1" / "S1.bam.bai").is_file())


if __name__ == "__main__":
    unittest.main()
