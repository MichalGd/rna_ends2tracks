from __future__ import annotations

import argparse
import json
import sys
from collections.abc import Callable
from pathlib import Path

from . import __version__
from .apa_b import apa_b
from .apa_mcell import active_pas_stage, apa_statistics_stage, exact_ends_stage
from .cleanup import clean_intermediates
from .compare import compare_apa
from .config import (
    ConfigError,
    RunPlan,
    build_conf_plan,
    sample_universe,
    workflow_requirements,
    write_plan,
)
from .dge import gene_expression
from .external import event
from .locking import run_lock
from .paths import workflow_asset
from .preprocess import preprocess
from .report import make_report
from .tracks import make_tracks

STEPS = (
    "validate", "alignment", "exact_ends", "active_pas", "gene_expression",
    "apa_a", "apa_b", "apa_comparison", "tracks", "report", "cleanup",
)
ALIASES = {"preprocess": "alignment", "dge": "gene_expression", "compare": "apa_comparison"}


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(
        prog="rna-ends2tracks",
        description="Config-driven QuantSeq REV gene-expression and APA workflow",
    )
    result.add_argument("--version", action="version", version=__version__)
    result.add_argument("--dry-run", action="store_true", help="Validate and record planned work without running tools")
    result.add_argument("--from-step", choices=(*STEPS, *ALIASES), help="Resume at this ordered step")
    result.add_argument("--stop-after", choices=(*STEPS, *ALIASES), help="Stop after this ordered step")
    result.add_argument("--force-step", action="append", default=[], choices=(*STEPS, *ALIASES),
                        help="Ignore matching receipt for this step; may be repeated")
    result.add_argument("--skip-input-checks", action="store_true",
                        help="Metadata-only validation; intended for portable CI examples")
    result.add_argument("config", help="Restricted KEY=value config.conf")
    return result


def _normal_step(value: str | None, fallback: str) -> str:
    return ALIASES.get(value or fallback, value or fallback)


def _selected_steps(args: argparse.Namespace) -> list[str]:
    first = STEPS.index(_normal_step(args.from_step, STEPS[0]))
    last = STEPS.index(_normal_step(args.stop_after, STEPS[-1]))
    if first > last:
        raise ConfigError("--from-step occurs after --stop-after")
    return list(STEPS[first:last + 1])


def _summary(plan: RunPlan, results: Path) -> dict[str, object]:
    genomes = sorted({sample["genome"] for sample in plan.samples})
    return {
        "status": "valid", "project_id": plan.project["project_id"], "genomes": genomes,
        "samples": len(plan.samples),
        "biological_replicates": len({(sample["genome"], sample["biological_replicate_id"])
                                      for sample in plan.samples}),
        "technical_replicates": len({(row["sample_id"], row["technical_replicate_id"])
                                     for row in plan.sample_rows}),
        "lanes": len(plan.sample_rows), "contrasts": len(plan.contrasts),
        "paired_contrasts": sum(row.get("design_mode") == "paired" for row in plan.contrasts),
        "unpaired_contrasts": sum(row.get("design_mode") == "unpaired" for row in plan.contrasts),
        "results": str(results),
    }


def execute(args: argparse.Namespace) -> int:
    if args.skip_input_checks and not args.dry_run and _normal_step(args.stop_after, "") != "validate":
        raise ConfigError("--skip-input-checks is allowed only with --dry-run or --stop-after validate")
    plan = build_conf_plan(args.config, check_inputs=not args.skip_input_checks)
    results = Path(plan.project["output_dir"]).expanduser().resolve()
    results.mkdir(parents=True, exist_ok=True)
    universe_path = results / "04_active_pas" / "sample_set.json"
    current_universe = sample_universe(plan)
    if universe_path.is_file():
        try:
            prior_universe = json.loads(universe_path.read_text(encoding="utf-8"))
        except (OSError, ValueError) as exc:
            raise ConfigError(f"Cannot validate existing PAS-universe signature: {universe_path}") from exc
        if prior_universe != current_universe:
            raise ConfigError(
                "OUTPUT_DIR already contains a different active-PAS sample universe; use a new OUTPUT_DIR"
            )
    steps = _selected_steps(args)
    forced = {_normal_step(value, value) for value in args.force_step}
    script_root = workflow_asset("scripts")

    with run_lock(results):
        metadata = results / "00_metadata"
        write_plan(plan, metadata)
        event(results / "logs", "validate", "completed",
              f"{len(plan.samples)} samples; {len(plan.contrasts)} within-genome contrasts")
        if steps == ["validate"]:
            print(json.dumps(_summary(plan, results), indent=2))
            return 0

        actions: dict[str, Callable[[], object]] = {
            "validate": lambda: None,
            "alignment": lambda: preprocess(plan, results, args.dry_run, "alignment" in forced),
            "exact_ends": lambda: exact_ends_stage(plan, results, args.dry_run, "exact_ends" in forced),
            "active_pas": lambda: active_pas_stage(plan, results, args.dry_run, "active_pas" in forced),
            "gene_expression": lambda: gene_expression(
                plan, results, script_root, args.dry_run, "gene_expression" in forced),
            "apa_a": lambda: apa_statistics_stage(plan, results, script_root, args.dry_run, "apa_a" in forced),
            "apa_b": lambda: apa_b(plan, results, script_root, args.dry_run, "apa_b" in forced),
            "apa_comparison": lambda: compare_apa(plan, results, force="apa_comparison" in forced),
            "tracks": lambda: make_tracks(plan, results, args.dry_run, "tracks" in forced),
            "report": lambda: make_report(plan, results, args.dry_run, "report" in forced),
            "cleanup": lambda: clean_intermediates(plan, results, args.dry_run, "cleanup" in forced),
        }
        modules = plan.project.get("modules", {})
        requirements = workflow_requirements(plan.project)
        for step in steps:
            if step == "exact_ends" and not requirements["exact_ends"]:
                event(results / "logs", step, "disabled", "No enabled module or track family requires exact ends")
                continue
            if step == "active_pas" and not requirements["active_pas"]:
                event(results / "logs", step, "disabled", "No enabled module or normalization requires active PAS")
                continue
            if step == "gene_expression" and not modules.get("gene_expression", True):
                event(results / "logs", step, "disabled", "RUN_GENE_EXPRESSION=false")
                continue
            if step == "apa_a" and not modules.get("apa_a", True):
                event(results / "logs", step, "disabled", "RUN_APA_A_MCELL2019=false")
                continue
            if step == "apa_b" and not plan.project["apa_b"]["enabled"]:
                event(results / "logs", step, "disabled", "RUN_APA_B=false")
                continue
            if step == "apa_comparison" and not requirements["apa_comparison"]:
                event(results / "logs", step, "disabled", "APA comparison requires both APA-A and APA-B")
                continue
            if step == "tracks" and not modules.get("tracks", True):
                event(results / "logs", step, "disabled", "RUN_TRACKS=false")
                continue
            # Reporting and destructive cleanup are meaningless during a dry run.
            if args.dry_run and step in {"report", "cleanup", "apa_comparison"}:
                event(results / "logs", step, "dry_run", "Would execute after successful upstream stages")
                continue
            actions[step]()
    print(json.dumps({**_summary(plan, results), "status": "dry_run" if args.dry_run else "completed",
                      "steps": steps}, indent=2))
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
