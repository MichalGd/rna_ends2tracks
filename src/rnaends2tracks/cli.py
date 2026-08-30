from __future__ import annotations

import argparse
import json
import os
import shutil
import sys
from collections.abc import Callable
from pathlib import Path

from . import __version__
from .apa_b import apa_b
from .apa_mcell import active_pas_stage, apa_statistics_stage, exact_ends_stage
from .cleanup import clean_intermediates
from .compare import compare_apa
from .conf import ConfError, project_from_conf
from .config import (
    ConfigError,
    RunPlan,
    build_conf_plan,
    sample_universe,
    workflow_requirements,
    write_plan,
)
from .dge import gene_expression
from .enrichment import enrichment
from .execution import run_bounded
from .external import event, read_run_status
from .locking import run_lock
from .paths import workflow_asset
from .preprocess import preprocess
from .report import make_report
from .rseqc import rseqc
from .tracks import make_c0_tracks, make_tracks

STEPS = (
    "validate", "alignment", "rseqc", "c0_tracks", "exact_ends", "active_pas", "gene_expression",
    "apa_a", "apa_b", "apa_comparison", "enrichment", "tracks", "report", "cleanup",
)
ALIASES = {
    "preprocess": "alignment", "early_tracks": "c0_tracks",
    "dge": "gene_expression", "compare": "apa_comparison",
}


def _downstream_branch_sequences(runnable: list[str]) -> list[tuple[str, tuple[str, ...]]]:
    """Return dependency-safe downstream branches in critical-path order."""
    branches: list[tuple[str, tuple[str, ...]]] = []
    remaining = list(runnable)
    if "gene_expression" in remaining and "tracks" in remaining:
        branches.append(("gene_expression_then_tracks", ("gene_expression", "tracks")))
        remaining.remove("gene_expression")
        remaining.remove("tracks")
    branches.extend((candidate, (candidate,)) for candidate in remaining)
    priority = {"apa_b": 0, "gene_expression_then_tracks": 1, "apa_a": 2}
    branches.sort(key=lambda item: priority.get(item[0], 3))
    return branches


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
    result.add_argument("--config", dest="config_option", metavar="FILE",
                        help="Restricted KEY=value config.conf")
    result.add_argument("config", nargs="?",
                        help="Restricted KEY=value config.conf (legacy positional form)")
    return result


def status_parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(
        prog="rna-ends2tracks status",
        description="Show the latest stage and status for a workflow configuration or result directory",
    )
    result.add_argument("target", help="config.conf or workflow result directory")
    result.add_argument("--json", action="store_true", help="Print the complete machine-readable status")
    return result


def show_status(target: str | Path, as_json: bool = False) -> int:
    path = Path(target).expanduser().resolve()
    if path.is_file():
        try:
            project, _samplesheet = project_from_conf(path)
        except ConfError as exc:
            raise ConfigError(str(exc)) from exc
        results = Path(project["output_dir"]).expanduser().resolve()
    else:
        results = path
    payload = read_run_status(results)
    observations = _status_observations(results, payload)
    if as_json:
        print(json.dumps({**payload, "observed": observations}, indent=2, sort_keys=True))
    else:
        print(f"Workflow status: {payload.get('workflow_status', 'unknown')}")
        print(f"Current stage:   {payload.get('current_stage', 'unknown')}")
        print(f"Last update:     {payload.get('updated_at', 'unknown')}")
        print(f"Message:         {payload.get('last_message', '')}")
        print(f"Workflow PID:    {observations['workflow_pid']} ({observations['process_state']})")
        print(f"Output directory:{' ' if str(results) else ''}{results}")
        print(f"Master log:      {results / 'rna_ends2tracks.log'}")
        print(f"Free disk:       {observations['free_disk_gb']:.1f} GiB")
        print("Outputs:         " + ", ".join(
            f"{key}={value}" for key, value in observations["outputs"].items()
        ))
        print("Stages:")
        stages = payload.get("stages", {})
        for step in STEPS:
            record = stages.get(step, {}) if isinstance(stages, dict) else {}
            reported = str(record.get("status", "pending"))
            receipt = observations["stage_receipts"].get(step, "missing")
            if reported in {"pending", "dry_run", "started"} and receipt == "completed":
                reported = "completed (receipt)"
            print(f"  {step:<18} {reported}")
    return 0


def _process_state(pid: object, workflow_status: object) -> str:
    if workflow_status in {"completed", "failed"}:
        return "not running"
    try:
        numeric_pid = int(pid)
        os.kill(numeric_pid, 0)
    except (TypeError, ValueError, ProcessLookupError):
        return "not running"
    except PermissionError:
        return "running (owned by another user)"
    except OSError:
        return "unknown"
    return "running"


def _table_rows(path: Path) -> int:
    if not path.is_file():
        return 0
    with path.open(encoding="utf-8-sig") as handle:
        return max(sum(1 for line in handle if line.strip()) - 1, 0)


def _receipt_completed(path: Path) -> bool:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return False
    return payload.get("schema_version") == 1 and payload.get("exit_status") == 0


def _observed_stage_receipts(results: Path) -> dict[str, str]:
    locations = {
        "alignment": results / "02_alignment" / "run_receipt.json",
        "rseqc": results / "01_qc" / "rseqc" / "run_receipt.json",
        "c0_tracks": results / "09_tracks" / ".stage_receipts" / "tracks_c0" / "run_receipt.json",
        "exact_ends": results / "03_exact_ends" / "run_receipt.json",
        "active_pas": results / "04_active_pas" / "run_receipt.json",
        "gene_expression": results / "05_gene_expression" / "run_receipt.json",
        "apa_a": results / "06_apa_a_mcell2019" / "run_receipt.json",
        "apa_b": results / "07_apa_b" / "run_receipt.json",
        "apa_comparison": results / "08_apa_comparison" / "run_receipt.json",
        "enrichment": results / "10_reports" / "enrichment_summary" / "run_receipt.json",
        "tracks": results / "09_tracks" / "run_receipt.json",
        "report": results / "10_reports" / "run_receipt.json",
        "cleanup": results / "provenance" / "cleanup" / "run_receipt.json",
    }
    return {
        step: "completed" if _receipt_completed(path) else "missing"
        for step, path in locations.items()
    }


def _status_observations(results: Path, payload: dict[str, object]) -> dict[str, object]:
    alignment = results / "02_alignment"
    final_bams = sum(
        (directory / f"{directory.name}.bam").is_file()
        for directory in alignment.iterdir()
        if alignment.is_dir() and directory.is_dir()
    ) if alignment.is_dir() else 0
    usage = shutil.disk_usage(results)
    pid = payload.get("workflow_pid", payload.get("pid", "unknown"))
    enrichment_index = results / "10_reports" / "enrichment_summary" / "enrichment_index.tsv"
    enrichment_jobs = _table_rows(enrichment_index)
    if not enrichment_jobs:
        enrichment_jobs = sum(
            _receipt_completed(path)
            for root in ("05_gene_expression", "06_apa_a_mcell2019", "07_apa_b")
            for path in (results / root).rglob("enrichment/*/run_receipt.json")
        )
    outputs = {
        "contrasts": _table_rows(results / "00_metadata" / "contrasts.tsv"),
        "BAMs": final_bams,
        "RSeQC": _table_rows(results / "01_qc" / "rseqc" / "rseqc_summary.tsv"),
        "BigWigs": len(list((results / "09_tracks").rglob("*.bw"))),
        "DGE": len(list((results / "05_gene_expression").rglob("*.deseq2.tsv"))),
        "APA-A": len(list((results / "06_apa_a_mcell2019").rglob("*.dexseq.tsv"))),
        "APA-B": len(list((results / "07_apa_b").rglob("*.drimseq_stager.tsv"))),
        "enrichment": enrichment_jobs,
        "reports": int((results / "10_reports" / "report.html").is_file()),
    }
    return {
        "workflow_pid": pid,
        "process_state": _process_state(pid, payload.get("workflow_status")),
        "free_disk_gb": usage.free / 1024**3,
        "outputs": outputs,
        "stage_receipts": _observed_stage_receipts(results),
    }


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


def _resolve_run_config(args: argparse.Namespace) -> str:
    positional = getattr(args, "config", None)
    option = getattr(args, "config_option", None)
    if positional and option:
        raise ConfigError("Specify the configuration once, using either --config FILE or positional FILE")
    value = option or positional
    if not value:
        raise ConfigError("A configuration is required; use --config FILE or positional FILE")
    return str(value)


def execute(args: argparse.Namespace) -> int:
    args.config = _resolve_run_config(args)
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
        event(results / "logs", "workflow", "started", f"Selected stages: {','.join(steps)}")
        event(results / "logs", "validate", "completed",
              f"{len(plan.samples)} samples; {len(plan.contrasts)} within-genome contrasts")
        if steps == ["validate"]:
            event(results / "logs", "workflow", "completed", "Validation-only run completed")
            print(json.dumps(_summary(plan, results), indent=2))
            return 0

        actions: dict[str, Callable[[], object]] = {
            "validate": lambda: None,
            "alignment": lambda: preprocess(plan, results, args.dry_run, "alignment" in forced),
            "rseqc": lambda: rseqc(plan, results, args.dry_run, "rseqc" in forced),
            "c0_tracks": lambda: make_c0_tracks(plan, results, args.dry_run, "c0_tracks" in forced),
            "exact_ends": lambda: exact_ends_stage(plan, results, args.dry_run, "exact_ends" in forced),
            "active_pas": lambda: active_pas_stage(plan, results, args.dry_run, "active_pas" in forced),
            "gene_expression": lambda: gene_expression(
                plan, results, script_root, args.dry_run, "gene_expression" in forced),
            "apa_a": lambda: apa_statistics_stage(plan, results, script_root, args.dry_run, "apa_a" in forced),
            "apa_b": lambda: apa_b(plan, results, script_root, args.dry_run, "apa_b" in forced),
            "apa_comparison": lambda: compare_apa(plan, results, force="apa_comparison" in forced),
            "enrichment": lambda: enrichment(
                plan, results, script_root, args.dry_run, "enrichment" in forced),
            "tracks": lambda: make_tracks(plan, results, args.dry_run, "tracks" in forced),
            "report": lambda: make_report(plan, results, args.dry_run, "report" in forced),
            "cleanup": lambda: clean_intermediates(plan, results, args.dry_run, "cleanup" in forced),
        }
        modules = plan.project.get("modules", {})
        requirements = workflow_requirements(plan.project)
        processed: set[str] = set()

        def disabled_reason(step: str) -> str | None:
            if step == "rseqc" and not modules.get("rseqc", False):
                return "RUN_RSEQC=false"
            if step == "c0_tracks" and not (
                modules.get("tracks", True)
                and plan.project["tracks"].get("early_c0", True)
                and plan.project["tracks"]["families"].get("all_reads", False)
            ):
                return "Requires RUN_TRACKS, GENERATE_EARLY_C0_TRACKS and GENERATE_ALL_READ_TRACKS"
            if step == "exact_ends" and not requirements["exact_ends"]:
                return "No enabled module or track family requires exact ends"
            if step == "active_pas" and not requirements["active_pas"]:
                return "No enabled module or normalization requires active PAS"
            if step == "gene_expression" and not modules.get("gene_expression", True):
                return "RUN_GENE_EXPRESSION=false"
            if step == "apa_a" and not modules.get("apa_a", True):
                return "RUN_APA_A_MCELL2019=false"
            if step == "apa_b" and not plan.project["apa_b"]["enabled"]:
                return "RUN_APA_B=false"
            if step == "apa_comparison" and not requirements["apa_comparison"]:
                return "APA comparison requires both APA-A and APA-B"
            if step == "enrichment" and not (
                modules.get("dge_enrichment", True) or modules.get("apa_enrichment", True)
            ):
                return "DGE and APA enrichment are disabled"
            if step == "tracks" and not modules.get("tracks", True):
                return "RUN_TRACKS=false"
            return None

        # Once active PAS are available, three independent branches can run:
        # DGE followed by final tracks, APA-A, and APA-B.  Tracks consume the
        # C4 size factors published by DGE but do not depend on either APA
        # statistical engine, so keeping them in the DGE branch avoids waiting
        # for the long APA-B critical path.
        downstream_group = ("gene_expression", "apa_a", "apa_b", "tracks")
        for step in steps:
            if step in processed:
                continue
            reason = disabled_reason(step)
            if reason:
                event(results / "logs", step, "disabled", reason)
                processed.add(step)
                continue
            if (
                not args.dry_run
                and step in downstream_group
                and plan.project["resources"]["downstream"]["parallel_modules"] > 1
            ):
                candidates = [
                    candidate for candidate in downstream_group
                    if candidate in steps and candidate not in processed
                ]
                runnable: list[str] = []
                for candidate in candidates:
                    reason = disabled_reason(candidate)
                    if reason:
                        event(results / "logs", candidate, "disabled", reason)
                    else:
                        runnable.append(candidate)
                    processed.add(candidate)
                branch_sequences = _downstream_branch_sequences(runnable)

                def module_worker(sequence: tuple[str, ...]) -> Callable[[], object]:
                    def worker() -> object:
                        value: object = None
                        for candidate in sequence:
                            event(results / "logs", candidate, "started", "Stage started")
                            try:
                                value = actions[candidate]()
                            except Exception as exc:
                                event(results / "logs", candidate, "failed", str(exc))
                                raise
                            event(results / "logs", candidate, "completed", "Stage completed")
                        return value
                    return worker

                if len(branch_sequences) > 1:
                    event(
                        results / "logs", "downstream", "started",
                        "Running dependency-aware downstream branches concurrently within the "
                        "aggregate resource ceiling: "
                        + ",".join(label for label, _sequence in branch_sequences),
                    )

                    try:
                        run_bounded(
                            "downstream_branches",
                            [
                                (label, module_worker(sequence))
                                for label, sequence in branch_sequences
                            ],
                            plan.project["resources"]["downstream"]["parallel_modules"],
                            results / ".checkpoints" / "timings" / "downstream_branches",
                        )
                    except Exception as exc:
                        event(results / "logs", "downstream", "failed", str(exc))
                        event(results / "logs", "workflow", "failed", f"Stopped during downstream branches: {exc}")
                        raise
                    event(results / "logs", "downstream", "completed", "Concurrent downstream modules completed")
                    continue
                if branch_sequences and len(branch_sequences[0][1]) > 1:
                    try:
                        module_worker(branch_sequences[0][1])()
                    except Exception as exc:
                        event(
                            results / "logs", "workflow", "failed",
                            f"Stopped during downstream branch {branch_sequences[0][0]}: {exc}",
                        )
                        raise
                    continue
                if runnable:
                    step = runnable[0]
            # Reporting and destructive cleanup are meaningless during a dry run.
            if args.dry_run and step in {"report", "cleanup", "apa_comparison"}:
                event(results / "logs", step, "dry_run", "Would execute after successful upstream stages")
                continue
            event(results / "logs", step, "started", "Stage started")
            try:
                actions[step]()
            except Exception as exc:
                event(results / "logs", step, "failed", str(exc))
                event(results / "logs", "workflow", "failed", f"Stopped during {step}: {exc}")
                raise
            event(results / "logs", step, "completed", "Stage completed")
            processed.add(step)
        event(results / "logs", "workflow", "completed",
              "Dry run completed" if args.dry_run else "All selected stages completed")
    print(json.dumps({**_summary(plan, results), "status": "dry_run" if args.dry_run else "completed",
                      "steps": steps}, indent=2))
    return 0


def main() -> None:
    try:
        if len(sys.argv) > 1 and sys.argv[1] == "status":
            status_args = status_parser().parse_args(sys.argv[2:])
            raise SystemExit(show_status(status_args.target, status_args.json))
        raise SystemExit(execute(parser().parse_args()))
    except (ConfigError, RuntimeError, OSError, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(2)


if __name__ == "__main__":
    main()
