import json
import os
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from rnaends2tracks.receipts import receipt_valid, write_receipt


class ReceiptTests(unittest.TestCase):
    def test_small_outputs_use_sha256_and_detect_same_size_changes(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary); output = root / "small.tsv"; output.write_bytes(b"abc")
            receipt_dir = root / "new" / "receipt"
            with patch("rnaends2tracks.receipts.HASH_LIMIT_BYTES", 100):
                receipt = write_receipt("test", receipt_dir, "sig", [output], ["test"])
                payload = json.loads(receipt.read_text(encoding="utf-8"))
                self.assertEqual(payload["outputs"][0]["validation"], "sha256")
                self.assertTrue(receipt_valid(receipt_dir, "sig"))
                output.write_bytes(b"abd")
                self.assertFalse(receipt_valid(receipt_dir, "sig"))

    def test_large_outputs_avoid_hashing_and_validate_size_mtime(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary); output = root / "large.bam"; output.write_bytes(b"abc")
            receipt_dir = root / "receipt"
            with patch("rnaends2tracks.receipts.HASH_LIMIT_BYTES", 1):
                receipt = write_receipt("test", receipt_dir, "sig", [output], ["test"])
                payload = json.loads(receipt.read_text(encoding="utf-8"))
                record = payload["outputs"][0]
                self.assertEqual(record["validation"], "size_mtime")
                self.assertNotIn("sha256", record)
                self.assertTrue(receipt_valid(receipt_dir, "sig"))
                stat = output.stat()
                os.utime(output, ns=(stat.st_atime_ns, stat.st_mtime_ns + 1_000_000))
                self.assertFalse(receipt_valid(receipt_dir, "sig"))

    def test_post2_accepts_only_audited_alpha9_family_receipts(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary); output = root / "result.tsv"; output.write_bytes(b"abc")
            receipt_dir = root / "receipt"
            receipt = write_receipt("test", receipt_dir, "sig", [output], ["test"])
            payload = json.loads(receipt.read_text(encoding="utf-8"))
            with patch("rnaends2tracks.receipts.__version__", "0.1.0a9.post2"):
                payload["workflow_version"] = "0.1.0a9"
                receipt.write_text(json.dumps(payload), encoding="utf-8")
                self.assertTrue(receipt_valid(receipt_dir, "sig"))
                payload["workflow_version"] = "0.1.0a9.post1"
                receipt.write_text(json.dumps(payload), encoding="utf-8")
                self.assertTrue(receipt_valid(receipt_dir, "sig"))
                payload["workflow_version"] = "0.1.0a8"
                receipt.write_text(json.dumps(payload), encoding="utf-8")
                self.assertFalse(receipt_valid(receipt_dir, "sig"))

    def test_alpha11_post1_accepts_only_alpha11_base_receipts(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary); output = root / "result.tsv"; output.write_bytes(b"abc")
            receipt_dir = root / "receipt"
            receipt = write_receipt("test", receipt_dir, "sig", [output], ["test"])
            payload = json.loads(receipt.read_text(encoding="utf-8"))
            with patch("rnaends2tracks.receipts.__version__", "0.1.0a11.post1"):
                payload["workflow_version"] = "0.1.0a11"
                receipt.write_text(json.dumps(payload), encoding="utf-8")
                self.assertTrue(receipt_valid(receipt_dir, "sig"))
                payload["workflow_version"] = "0.1.0a10.post7"
                receipt.write_text(json.dumps(payload), encoding="utf-8")
                self.assertFalse(receipt_valid(receipt_dir, "sig"))

    def test_alpha11_post2_accepts_audited_alpha11_receipts(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary); output = root / "result.tsv"; output.write_bytes(b"abc")
            receipt_dir = root / "receipt"
            receipt = write_receipt("test", receipt_dir, "sig", [output], ["test"])
            payload = json.loads(receipt.read_text(encoding="utf-8"))
            with patch("rnaends2tracks.receipts.__version__", "0.1.0a11.post2"):
                for version in ("0.1.0a11", "0.1.0a11.post1"):
                    payload["workflow_version"] = version
                    receipt.write_text(json.dumps(payload), encoding="utf-8")
                    self.assertTrue(receipt_valid(receipt_dir, "sig"))
                payload["workflow_version"] = "0.1.0a10.post7"
                receipt.write_text(json.dumps(payload), encoding="utf-8")
                self.assertFalse(receipt_valid(receipt_dir, "sig"))

    def test_alpha11_post3_accepts_audited_alpha11_receipts(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary); output = root / "result.tsv"; output.write_bytes(b"abc")
            receipt_dir = root / "receipt"
            receipt = write_receipt("test", receipt_dir, "sig", [output], ["test"])
            payload = json.loads(receipt.read_text(encoding="utf-8"))
            with patch("rnaends2tracks.receipts.__version__", "0.1.0a11.post3"):
                for version in ("0.1.0a11", "0.1.0a11.post1", "0.1.0a11.post2"):
                    payload["workflow_version"] = version
                    receipt.write_text(json.dumps(payload), encoding="utf-8")
                    self.assertTrue(receipt_valid(receipt_dir, "sig"))
                payload["workflow_version"] = "0.1.0a10.post7"
                receipt.write_text(json.dumps(payload), encoding="utf-8")
                self.assertFalse(receipt_valid(receipt_dir, "sig"))

    def test_alpha12_accepts_audited_alpha11_receipts(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary); output = root / "result.tsv"; output.write_bytes(b"abc")
            receipt_dir = root / "receipt"
            receipt = write_receipt("test", receipt_dir, "sig", [output], ["test"])
            payload = json.loads(receipt.read_text(encoding="utf-8"))
            with patch("rnaends2tracks.receipts.__version__", "0.1.0a12"):
                for version in (
                    "0.1.0a11", "0.1.0a11.post1", "0.1.0a11.post2", "0.1.0a11.post3",
                ):
                    payload["workflow_version"] = version
                    receipt.write_text(json.dumps(payload), encoding="utf-8")
                    self.assertTrue(receipt_valid(receipt_dir, "sig"))
                payload["workflow_version"] = "0.1.0a10.post7"
                receipt.write_text(json.dumps(payload), encoding="utf-8")
                self.assertFalse(receipt_valid(receipt_dir, "sig"))


if __name__ == "__main__":
    unittest.main()
