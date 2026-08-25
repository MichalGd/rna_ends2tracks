from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from . import __version__
from .apa_a import apa_a
from .apa_b import apa_b
from .compare import compare_apa
from .config import ConfigError, build_plan, write_plan
from .dge import gene_expression
from .external import event
from .preprocess import preprocess
from .paths import workflow_asset
from .report import make_report
from .tracks import make_tracks


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(prog="rna-ends2tracks", description="QuantSeq REV DGE and independent APA analyses")
    result.add_argument("--version", action="version", version=__version__)
    result.add_argument("--config", required=True, help="Project YAML")
    result.add_argument("--samplesheet", required=True, help="Lane-level CSV samplesheet")
    result.add_argument("--results", help="Results directory; overrides output_dir in project YAML")
    result.add_argument("--skip-input-checks", action="store_true", help="Validate metadata without requiring FASTQ/reference files")
    result.add_argument("--dry-run", action="store_true", help="Write the plan and log commands without running tools")
    result.add_argument("--force-module", action="store_true", help="Ignore a matching module receipt")
    sub = result.add_subparsers(dest="module", required=True)
    for name in ("validate", "preprocess", "dge", "apa-a", "apa-b", "compare", "tracks", "report", "all"):
        child = sub.add_parser(name)
        if name == "compare":
            child.add_argument("--tolerance", type=int, default=24)
    return result


def _results_dir(args: argparse.Namespace, project: dict) -> Path:
    if args.results:
        return Path(args.results).expanduser().resolve()
    configured = project.get("output_dir", "results")
    path = Path(str(configured)).expanduser()
    if not path.is_absolute():
        path = Path(project["_config_path"]).parent / path
    return path.resolve()


def execute(args: argparse.Namespace) -> int:
    check_inputs = not args.skip_input_checks
    plan = build_plan(args.config, args.samplesheet, check_inputs=check_inputs)
    results = _results_dir(args, plan.project)
    results.mkdir(parents=True, exist_ok=True)
    metadata = results / "00_metadata"
    write_plan(plan, metadata)
    event(results / "provenance" / "logs", "validate", "completed", f"{len(plan.samples)} samples; {len(plan.contrasts)} contrasts")
    script_root = workflow_asset("scripts")
    if args.module == "validate":
        print(json.dumps({
            "status": "valid", "species": plan.reference["species"], "assembly": plan.reference["assembly"],
            "samples": len(plan.samples), "lanes": len(plan.sample_rows), "contrasts": len(plan.contrasts),
            "results": str(results),
            "warnings": plan.reference.get("_warnings", []),
        }, indent=2))
        return 0
    if args.module in {"preprocess", "all"}:
        preprocess(plan, results, args.dry_run, args.force_module)
    if args.module in {"dge", "all"}:
        gene_expression(plan, results, script_root, args.dry_run, args.force_module)
    if args.module in {"apa-a", "all"}:
        apa_a(plan, results, script_root, args.dry_run, args.force_module)
    if args.module in {"apa-b", "all"}:
        apa_b(plan, results, script_root, args.dry_run, args.force_module)
    if args.module in {"tracks", "all"}:
        make_tracks(plan, results, args.dry_run, args.force_module)
    if args.module == "compare" or (args.module == "all" and plan.project.get("apa_b", {}).get("enabled", False)):
        if not args.dry_run:
            compare_apa(results, getattr(args, "tolerance", 24), args.force_module)
    if args.module in {"report", "all"} and not args.dry_run:
        print(make_report(results, args.force_module))
    return 0


def main() -> None:
    args = parser().parse_args()
    try:
        raise SystemExit(execute(args))
    except (ConfigError, RuntimeError, OSError, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(2)


if __name__ == "__main__":
    main()
