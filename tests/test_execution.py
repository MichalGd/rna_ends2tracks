import os
import tempfile
import unittest
from pathlib import Path

from rnaends2tracks.execution import run_bounded_processes


def process_identity(value):
    return os.getpid(), value


def process_failure():
    raise ValueError("intentional worker failure")


class ProcessExecutionTests(unittest.TestCase):
    def test_process_pool_preserves_order_and_records_worker_pids(self):
        with tempfile.TemporaryDirectory() as temporary:
            timing = Path(temporary)
            jobs = [
                ("second", process_identity, (2,)),
                ("first", process_identity, (1,)),
            ]
            results = run_bounded_processes("test_processes", jobs, 2, timing)
            self.assertEqual([value for _pid, value in results], [2, 1])
            self.assertTrue(all(pid != os.getpid() for pid, _value in results))
            self.assertIn('"executor": "process"', (timing / "second.json").read_text())

    def test_process_pool_reports_labelled_failure(self):
        with tempfile.TemporaryDirectory() as temporary:
            with self.assertRaisesRegex(RuntimeError, "broken: intentional worker failure"):
                run_bounded_processes(
                    "test_processes", [("broken", process_failure, ())], 1, Path(temporary)
                )


if __name__ == "__main__":
    unittest.main()
