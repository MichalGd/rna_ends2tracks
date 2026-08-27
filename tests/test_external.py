import json
import sys
import tempfile
import unittest
from contextlib import redirect_stdout
from io import StringIO
from pathlib import Path

from rnaends2tracks.cli import show_status
from rnaends2tracks.external import event, run


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

    def test_status_command_accepts_results_directory(self):
        with tempfile.TemporaryDirectory() as temporary:
            results = Path(temporary) / "results"
            event(results / "logs", "workflow", "completed", "Done")
            output = StringIO()
            with redirect_stdout(output):
                self.assertEqual(show_status(results), 0)
            self.assertIn("Workflow status: completed", output.getvalue())
            self.assertIn(str(results / "rna_ends2tracks.log"), output.getvalue())


if __name__ == "__main__":
    unittest.main()
