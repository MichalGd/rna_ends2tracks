import os
import tempfile
import unittest
from functools import partial
from pathlib import Path

from rnaends2tracks.execution import DEFAULT_RESOURCES, resource_plan_rows, run_bounded_processes


def process_identity(value):
    return os.getpid(), value


def process_failure():
    raise ValueError("intentional worker failure")


def append_value(target, value):
    target.append(value)
    return value


class ProcessExecutionTests(unittest.TestCase):
    def test_resource_plan_declares_exact_end_process_executor(self):
        rows = resource_plan_rows(DEFAULT_RESOURCES, {"samples": 4})
        exact = next(row for row in rows if row["work_unit"] == "exact_end_extraction")
        self.assertEqual(exact["executor"], "python_process")
        tracks = [row for row in rows if row["stage"] == "tracks"]
        self.assertEqual({row["work_unit"] for row in tracks}, {"c0_sample", "end_sample"})
        self.assertTrue(all(row["executor"] == "python_process_and_external_process" for row in tracks))
        apa_b = [row for row in rows if row["stage"] == "apa_b"]
        self.assertEqual(
            {row["work_unit"] for row in apa_b},
            {"endpoint_preparation", "polyaseqtrap_cluster", "deepip", "contrast"},
        )
        downstream = next(row for row in rows if row["work_unit"] == "module_overlap")
        self.assertEqual(downstream["executor"], "thread_coordinator_with_nested_bounded_pools")
        self.assertEqual(downstream["units"], 3)
        self.assertEqual(downstream["max_threads"], 8)
        self.assertEqual(downstream["max_memory_gb"], 32)
        self.assertEqual(downstream["budget_status"], "PASS")

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

    def test_thread_pool_reports_each_terminal_worker_state(self):
        from rnaends2tracks.execution import run_bounded

        with tempfile.TemporaryDirectory() as temporary:
            returned = []
            progress = []
            completed = []
            jobs = [
                ("one", partial(append_value, returned, 1)),
                ("two", partial(append_value, returned, 2)),
            ]
            self.assertEqual(
                run_bounded(
                    "threads", jobs, 2, Path(temporary),
                    progress=lambda label, status: progress.append((label, status)),
                    on_completed=lambda label, value: completed.append((label, value)),
                ),
                [1, 2],
            )
            self.assertEqual(set(progress), {("one", "completed"), ("two", "completed")})
            self.assertEqual(set(completed), {("one", 1), ("two", 2)})


if __name__ == "__main__":
    unittest.main()
