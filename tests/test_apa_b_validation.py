import json
import tempfile
import unittest
from pathlib import Path

from rnaends2tracks.apa_b import _validate_engine_provenance, _validation_manifest


def accepted_manifest():
    return {
        "schema_version": 1, "status": "accepted",
        "workflow_adapter": {"release": "v0.1.0-alpha.10.post4", "source_commit": "7654321"},
        "engine": {"name": "PolyAseqTrap adapter", "source_commit": "1234567"},
        "model": {"sha256": "a" * 64}, "environment": {"sha256": "b" * 64},
        "umi_present": False, "coordinate_deduplication": False,
        "quantseq_rev_adaptation": "genomewide_no_tail_weighted_PAC",
        "library_protocols": ["quantseq_rev_v2_se"], "assemblies": ["GRCh38", "GRCm39"],
        "pilot": {"synthetic_pass": True, "c1_c1s_reuse_equivalent": True,
                  "real_quantseq_rev_canaries": {"GRCh38": "PASS", "GRCm39": "PASS"}},
        "reviewed_by": "reviewer", "accepted_at": "2026-08-27T00:00:00Z",
    }


class ApaBValidationTests(unittest.TestCase):
    def test_accepted_manifest_and_matching_engine_provenance(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary); manifest = root / "validation.json"; provenance = root / "engine.json"
            value = accepted_manifest(); manifest.write_text(json.dumps(value), encoding="utf-8")
            observed = {
                "assembly": "GRCm39", "workflow_adapter": value["workflow_adapter"],
                "engine": value["engine"], "model": value["model"],
                "environment": value["environment"], "umi_present": False,
                "coordinate_deduplication": False, "library_protocol": "quantseq_rev_v2_se",
            }
            provenance.write_text(json.dumps(observed), encoding="utf-8")
            accepted = _validation_manifest(manifest, ["GRCm39"])
            _validate_engine_provenance(provenance, accepted, "GRCm39")

    def test_paired_protocol_requires_explicit_paired_acceptance(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "validation.json"; value = accepted_manifest()
            path.write_text(json.dumps(value), encoding="utf-8")
            with self.assertRaisesRegex(RuntimeError, "quantseq_rev_v2_pe"):
                _validation_manifest(path, ["GRCm39"], "quantseq_rev_v2_pe")
            value["library_protocols"].append("quantseq_rev_v2_pe")
            path.write_text(json.dumps(value), encoding="utf-8")
            self.assertEqual(
                _validation_manifest(path, ["GRCm39"], "quantseq_rev_v2_pe")["status"], "accepted"
            )

    def test_legacy_manifest_is_compatible_only_with_original_v2_se_protocol(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "validation.json"; value = accepted_manifest()
            value.pop("library_protocols")
            path.write_text(json.dumps(value), encoding="utf-8")
            observed = _validation_manifest(path, ["GRCm39"], "quantseq_rev_v2_se")
            self.assertEqual(observed["library_protocols"], ["quantseq_rev_v2_se"])
            with self.assertRaisesRegex(RuntimeError, "quantseq_rev_v2_pe"):
                _validation_manifest(path, ["GRCm39"], "quantseq_rev_v2_pe")

    def test_draft_or_missing_real_canary_is_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "validation.json"; value = accepted_manifest()
            value["pilot"]["real_quantseq_rev_canaries"]["GRCm39"] = "NOT_RUN"
            path.write_text(json.dumps(value), encoding="utf-8")
            with self.assertRaisesRegex(RuntimeError, "real QuantSeq REV canary"):
                _validation_manifest(path, ["GRCm39"])

    def test_pre_post4_manifest_without_reuse_equivalence_is_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "validation.json"; value = accepted_manifest()
            value["pilot"].pop("c1_c1s_reuse_equivalent")
            path.write_text(json.dumps(value), encoding="utf-8")
            with self.assertRaisesRegex(RuntimeError, "C1.C1S equivalence"):
                _validation_manifest(path, ["GRCm39"])


if __name__ == "__main__":
    unittest.main()
