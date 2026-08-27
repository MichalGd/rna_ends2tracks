import json
import sys
import tempfile
import unittest
from contextlib import redirect_stdout
from io import StringIO
from pathlib import Path

from rnaends2tracks.cli import show_status
from rnaends2tracks.external import event, progress_events, run, run_capture, run_to_path


def process_event(results, index):
    event(Path(results) / "logs", "parallel_test", "completed", f"worker={index}")
    return index


class UnifiedLoggingTests(unittest.TestCase):
    def test_event_updates_jsonl_master_log_and_status(self):
        with tempfile.TemporaryDirectory() as temporary:
            results = Path(temporary) / "results"
            event(results / "logs", "alignment", "started", "Aligning two lanes")
            self.assertIn('"module": "alignment"', (results / "logs" / "events.jsonl").read_text())
            self.assertIn("[alignment] STARTED Aligning two lanes", (results / "rna_ends2tracks.log").read_text())
            status = json.loads((results / "00_metadata" / "run_status.json").read_text())
            self.assertEqual(status["current_stage"], "alignment")
            self.assertEqual(status["workflow_status"], "running")

    def test_external_command_lifecycle_is_visible_in_master_log(self):
        with tempfile.TemporaryDirectory() as temporary:
            results = Path(temporary) / "results"
            detail = results / "logs" / "alignment" / "example.log"
            run([sys.executable, "-c", "print('ok')"], detail)
            master = (results / "rna_ends2tracks.log").read_text()
            self.assertIn("[alignment] STARTED", master)
            self.assertIn("[alignment] COMPLETED example: exit_status=0", master)
            self.assertIn("ok", detail.read_text())

    def test_capture_and_stdout_file_commands_are_logged(self):
        with tempfile.TemporaryDirectory() as temporary:
            results = Path(temporary) / "results"
            capture_log = results / "logs" / "test" / "capture.log"
            self.assertEqual(
                run_capture([sys.executable, "-c", "print('captured')"], capture_log).strip(),
                "captured",
            )
            output = results / "data.txt"
            run_to_path([sys.executable, "-c", "print('written')"], output,
                        results / "logs" / "test" / "write.log")
            self.assertEqual(output.read_text().strip(), "written")
            master = (results / "rna_ends2tracks.log").read_text()
            self.assertIn("capture: exit_status=0", master)
            self.assertIn("write: exit_status=0", master)

    def test_status_command_accepts_results_directory(self):
        with tempfile.TemporaryDirectory() as temporary:
            results = Path(temporary) / "results"
            metadata = results / "00_metadata"
            metadata.mkdir(parents=True)
            (metadata / "contrasts.tsv").write_text(
                "contrast_id\tnumerator\tdenominator\nB_vs_A\tB\tA\n", encoding="utf-8"
            )
            sample_dir = results / "02_alignment" / "S1"
            sample_dir.mkdir(parents=True)
            (sample_dir / "S1.bam").write_bytes(b"bam")
            track_dir = results / "09_tracks"
            track_dir.mkdir(parents=True)
            (track_dir / "S1.plus.bw").write_bytes(b"bw")
            event(results / "logs", "workflow", "completed", "Done")
            output = StringIO()
            with redirect_stdout(output):
                self.assertEqual(show_status(results), 0)
            self.assertIn("Workflow status: completed", output.getvalue())
            self.assertIn(str((results / "rna_ends2tracks.log").resolve()), output.getvalue())
            self.assertIn("contrasts=1, BAMs=1, BigWigs=1", output.getvalue())
            self.assertIn("Workflow PID:", output.getvalue())

    def test_status_json_adds_live_observations(self):
        with tempfile.TemporaryDirectory() as temporary:
            results = Path(temporary) / "results"
            event(results / "logs", "workflow", "started", "Running")
            output = StringIO()
            with redirect_stdout(output):
                self.assertEqual(show_status(results, as_json=True), 0)
            payload = json.loads(output.getvalue())
            self.assertEqual(payload["observed"]["process_state"], "running")
            self.assertEqual(payload["observed"]["workflow_pid"], payload["workflow_pid"])

    def test_event_updates_are_process_safe(self):
        from rnaends2tracks.execution import run_bounded_processes

        with tempfile.TemporaryDirectory() as temporary:
            results = Path(temporary) / "results"
            jobs = [(f"worker_{index}", process_event, (str(results), index)) for index in range(4)]
            self.assertEqual(
                run_bounded_processes("parallel_log", jobs, 4, results / "timings"),
                list(range(4)),
            )
            lines = (results / "logs" / "events.jsonl").read_text().splitlines()
            self.assertEqual(len(lines), 4)
            self.assertTrue(all(json.loads(line)["module"] == "parallel_test" for line in lines))
            json.loads((results / "00_metadata" / "run_status.json").read_text())

    def test_progress_callback_reports_completed_total_and_eta(self):
        with tempfile.TemporaryDirectory() as temporary:
            results = Path(temporary) / "results"
            progress = progress_events(results / "logs", "apa_a", 2, "contrast")
            progress("B_vs_A", "completed")
            progress("C_vs_A", "completed")
            master = (results / "rna_ends2tracks.log").read_text(encoding="utf-8")
            self.assertIn("contrast B_vs_A completed (1/2)", master)
            self.assertIn("contrast C_vs_A completed (2/2)", master)
            self.assertIn("ETA~00:00:00", master)


if __name__ == "__main__":
    unittest.main()
