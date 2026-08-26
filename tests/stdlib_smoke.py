"""Dependency-free contract smoke test for constrained build hosts.

Production and full tests use PyYAML and pysam from environment.yml. This file stubs only
their imports so core metadata mathematics and coordinate conventions can be checked with
a bare CPython runtime.
"""
from __future__ import annotations

import json
import sys
import types
from pathlib import Path
from typing import ClassVar

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

yaml_stub = types.ModuleType("yaml")
yaml_stub.safe_load = lambda handle: json.load(handle)
sys.modules.setdefault("yaml", yaml_stub)
pysam_stub = types.ModuleType("pysam")
pysam_stub.AlignedSegment = object
sys.modules.setdefault("pysam", pysam_stub)

from rnaends2tracks.apa_a import reverse_complement
from rnaends2tracks.config import (
    ConfigError,
    generate_contrasts,
    resolve_contrast_designs,
    validate_design,
)
from rnaends2tracks.mcell2019 import transcript_end


class Read:
    is_reverse = True
    reference_start = 100
    reference_end = 150
    cigartuples: ClassVar[list[tuple[int, int]]] = [(0, 50)]


def sample(identifier: str, condition: str, batch: str) -> dict[str, str]:
    return {
        "sample_id": identifier, "biological_replicate_id": identifier,
        "condition": condition, "batch": batch, "subject": "",
    }


assert transcript_end(Read()) == (149, "+", False)
assert reverse_complement("AACGT") == "ACGTT"
balanced = [sample("A1", "A", "X"), sample("A2", "A", "Y"), sample("B1", "B", "X"), sample("B2", "B", "Y")]
validate_design(balanced, "~ batch + condition")
assert generate_contrasts(balanced, ["A", "B"])[0]["contrast_id"] == "B_vs_A"
paired = [sample("A1", "A", "X"), sample("A2", "A", "X"), sample("B1", "B", "X"), sample("B2", "B", "X")]
paired[0]["subject"] = paired[2]["subject"] = "S1"
paired[1]["subject"] = paired[3]["subject"] = "S2"
resolved = resolve_contrast_designs(
    paired, generate_contrasts(paired, ["A", "B"]),
    {"design": "~ condition", "statistics": {"pairing": {"mode": "auto"}}},
)[0]
assert resolved["design_mode"] == "paired"
assert resolved["resolved_design"] == "~ subject + condition"
confounded = [sample("A1", "A", "X"), sample("A2", "A", "X"), sample("B1", "B", "Y"), sample("B2", "B", "Y")]
try:
    validate_design(confounded, "~ batch + condition")
except ConfigError:
    pass
else:
    raise AssertionError("Confounded design was accepted")
print("stdlib smoke: PASS")
