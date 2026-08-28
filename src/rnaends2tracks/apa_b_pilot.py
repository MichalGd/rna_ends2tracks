from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path

from .polyaseqtrap_adapter import _load_json


def _accepted_synthetic(path: Path) -> None:
    audit = _load_json(path, "APA-B synthetic pilot audit")
    required = {
        "status": "PASS",
        "coordinate_and_strand": True,
        "record_count_conserved": True,
        "duplicate_flagged_records_retained": True,
        "deepip_artifact_rejected": True,
        "deepip_true_pas_retained": True,
        "intragenic_site_retained": True,
    }
    failures = [key for key, expected in required.items() if audit.get(key) != expected]
    if failures:
        raise RuntimeError("Synthetic APA-B pilot has not passed: " + ", ".join(failures))


def _accepted_real_canary(assembly: str, directory: Path) -> None:
    required = ["pas_catalog.tsv", "pas_counts.tsv", "deepip_audit.tsv", "engine_provenance.json", "adapter_audit.json"]
    missing = [name for name in required if not (directory / name).is_file()]
    if missing:
        raise RuntimeError(f"{assembly} real canary is incomplete: " + ", ".join(missing))
    audit = _load_json(directory / "adapter_audit.json", f"{assembly} adapter audit")
    if int(audit.get("final_pas", 0)) < 1:
        raise RuntimeError(f"{assembly} real canary produced no retained PAS")
    for sample, record in audit.get("samples", {}).items():
        if record.get("eligible_records") != record.get("records_written"):
            raise RuntimeError(f"{assembly} real canary changed record count for {sample}")
    provenance = _load_json(directory / "engine_provenance.json", f"{assembly} engine provenance")
    if provenance.get("assembly") != assembly or provenance.get("coordinate_deduplication") is not False:
        raise RuntimeError(f"{assembly} real canary provenance violates the adapter contract")


def execute(args: argparse.Namespace) -> int:
    installation = _load_json(Path(args.installation_manifest), "APA-B installation manifest")
    _accepted_synthetic(Path(args.synthetic_audit))
    canaries: dict[str, str] = {}
    for specification in args.real_canary:
        assembly, separator, directory = specification.partition("=")
        if not separator or assembly not in {"GRCh38", "GRCm39"}:
            raise RuntimeError("--real-canary must use GRCh38=/path or GRCm39=/path")
        _accepted_real_canary(assembly, Path(directory))
        canaries[assembly] = "PASS"
    if not canaries:
        raise RuntimeError("At least one real QuantSeq REV canary is required")
    payload = {
        "schema_version": 1,
        "status": "accepted",
        "engine": installation["engine"],
        "deepip": {"source_commit": installation["deepip"]["source_commit"]},
        "models": {species: {"name": "DeepIP", "sha256": record["sha256"]}
                   for species, record in installation["models"].items()},
        # Compatibility for the main workflow's schema-v1 single-model validator.
        "model": {"name": "DeepIP", "sha256": next(iter(installation["models"].values()))["sha256"]},
        "environment": {"sha256": installation["environment"]["sha256"]},
        "assemblies": sorted(canaries),
        "library_protocols": ["quantseq_rev_v2_se"],
        "umi_present": False,
        "coordinate_deduplication": False,
        "quantseq_rev_adaptation": "genomewide_no_tail_weighted_PAC",
        "pilot": {"synthetic_pass": True, "real_quantseq_rev_canaries": canaries},
        "reviewed_by": args.reviewed_by,
        "accepted_at": datetime.now(timezone.utc).isoformat(),
    }
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return 0


def main() -> None:
    parser = argparse.ArgumentParser(description="Create an APA-B acceptance manifest from completed pilot audits")
    parser.add_argument("--installation-manifest", required=True)
    parser.add_argument("--synthetic-audit", required=True)
    parser.add_argument("--real-canary", action="append", default=[])
    parser.add_argument("--reviewed-by", required=True)
    parser.add_argument("--output", required=True)
    raise SystemExit(execute(parser.parse_args()))


if __name__ == "__main__":
    main()
