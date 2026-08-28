from pathlib import Path
import unittest

import yaml


ROOT = Path(__file__).resolve().parents[1]


class ApaBEnvironmentTests(unittest.TestCase):
    def test_tensorflow_is_pip_cpu_runtime_separate_from_conda_r_stack(self):
        payload = yaml.safe_load((ROOT / "environment.apa_b.yml").read_text(encoding="utf-8"))
        dependencies = payload["dependencies"]
        scalar = [entry for entry in dependencies if isinstance(entry, str)]
        pip_blocks = [entry["pip"] for entry in dependencies
                      if isinstance(entry, dict) and "pip" in entry]

        self.assertIn("python=3.10", scalar)
        self.assertIn("r-base=4.3", scalar)
        self.assertFalse(any(entry.startswith(("tensorflow=", "keras=")) for entry in scalar))
        self.assertEqual(len(pip_blocks), 1)
        self.assertIn("tensorflow-cpu==2.10.1", pip_blocks[0])
        self.assertIn("keras==2.10.0", pip_blocks[0])


if __name__ == "__main__":
    unittest.main()
