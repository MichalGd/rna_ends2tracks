from __future__ import annotations

import csv
import json
import os
import time
from collections.abc import Callable
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Any, TypeVar

T = TypeVar("T")
ProcessJob = tuple[str, Callable[..., T], tuple[Any, ...]]


DEFAULT_RESOURCES: dict[str, Any] = {
    "total_threads": 8,
    "total_memory_gb": 32,
    "temporary_directory": "",
    "preprocess": {
        "trim_parallel_jobs": 1,
        "star_parallel_jobs": 1,
        "fastqc_threads": 2,
        "bbduk_threads": 4,
        "bbduk_memory_gb": 8,
        "star_threads": 8,
        "star_memory_gb": 24,
        "samtools_threads": 4,
        "samtools_sort_memory_per_thread_gb": 2,
        "merge_parallel_jobs": 1,
        "merge_memory_gb": 8,
    },
    "dge": {
        "featurecounts_threads": 8,
        "featurecounts_memory_gb": 16,
        "contrast_parallel_jobs": 1,
        "contrast_threads": 1,
        "contrast_memory_gb": 16,
    },
    "apa_a": {
        "extraction_parallel_jobs": 2,
        "extraction_threads": 1,
        "extraction_memory_gb": 4,
        "contrast_parallel_jobs": 1,
        "contrast_threads": 1,
        "contrast_memory_gb": 16,
    },
    "apa_b": {
        "engine_threads": 8,
        "engine_memory_gb": 24,
        "contrast_parallel_jobs": 1,
        "contrast_threads": 1,
        "contrast_memory_gb": 16,
    },
    "enrichment": {
        "parallel_jobs": 1,
        "threads": 1,
        "memory_gb": 16,
    },
    "tracks": {
        "parallel_jobs": 2,
        "samtools_threads": 4,
        "memory_gb": 8,
    },
}


def _positive_int(value: Any, name: str) -> int:
    if isinstance(value, bool):
        raise ValueError(f"{name} must be a positive integer")  # noqa: TRY004
    try:
        result = int(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be a positive integer") from exc
    if result < 1 or str(value).strip() not in {str(result), f"{result}.0"}:
        raise ValueError(f"{name} must be a positive integer")
    return result


def _copy_defaults() -> dict[str, Any]:
    return json.loads(json.dumps(DEFAULT_RESOURCES))


def resolve_resources(project: dict[str, Any]) -> tuple[dict[str, Any], list[dict[str, str]]]:
    """Return a complete resource configuration and validate aggregate budgets."""
    supplied = project.get("resources", {})
    if not isinstance(supplied, dict):
        raise ValueError("resources must be a YAML mapping")  # noqa: TRY004
    resolved = _copy_defaults()
    warnings: list[dict[str, str]] = []

    legacy_threads = supplied.get("threads")
    legacy_memory = supplied.get("memory_gb")
    if legacy_threads is not None or legacy_memory is not None:
        warnings.append({
            "warning_code": "LEGACY_RESOURCE_CONFIGURATION",
            "message": (
                "resources.threads/resources.memory_gb were translated to alpha.6 global ceilings; "
                "migrate to total_threads, total_memory_gb, and stage-specific settings"
            ),
        })
        if "total_threads" not in supplied and legacy_threads is not None:
            resolved["total_threads"] = legacy_threads
        if "total_memory_gb" not in supplied and legacy_memory is not None:
            resolved["total_memory_gb"] = legacy_memory

    allowed_top = set(resolved) | {"threads", "memory_gb"}
    unknown_top = sorted(set(supplied) - allowed_top)
    if unknown_top:
        raise ValueError("Unknown resources settings: " + ", ".join(unknown_top))
    for key in ("total_threads", "total_memory_gb", "temporary_directory"):
        if key in supplied:
            resolved[key] = supplied[key]
    for stage in ("preprocess", "dge", "apa_a", "apa_b", "enrichment", "tracks"):
        stage_value = supplied.get(stage, {})
        if not isinstance(stage_value, dict):
            raise ValueError(f"resources.{stage} must be a YAML mapping")  # noqa: TRY004
        unknown = sorted(set(stage_value) - set(resolved[stage]))
        if unknown:
            raise ValueError(f"Unknown resources.{stage} settings: " + ", ".join(unknown))
        resolved[stage].update(stage_value)

    resolved["total_threads"] = _positive_int(resolved["total_threads"], "resources.total_threads")
    resolved["total_memory_gb"] = _positive_int(resolved["total_memory_gb"], "resources.total_memory_gb")
    resolved["temporary_directory"] = str(resolved.get("temporary_directory", "") or "")
    for stage in ("preprocess", "dge", "apa_a", "apa_b", "enrichment", "tracks"):
        for key, value in list(resolved[stage].items()):
            resolved[stage][key] = _positive_int(value, f"resources.{stage}.{key}")

    violations = [row for row in resource_plan_rows(resolved) if row["budget_status"] != "PASS"]
    if violations:
        details = "; ".join(
            f"{row['stage']}/{row['work_unit']}: {row['max_threads']} threads, "
            f"{row['max_memory_gb']} GB"
            for row in violations
        )
        raise ValueError(
            f"Resource plan exceeds global ceiling ({resolved['total_threads']} threads, "
            f"{resolved['total_memory_gb']} GB): {details}"
        )
    return resolved, warnings


def resource_plan_rows(resources: dict[str, Any], counts: dict[str, int] | None = None) -> list[dict[str, Any]]:
    counts = counts or {}
    pre = resources["preprocess"]
    dge = resources["dge"]
    apa_a = resources["apa_a"]
    apa_b = resources["apa_b"]
    enrichment = resources["enrichment"]
    tracks = resources["tracks"]
    definitions = [
        ("preprocess", "qc_and_trim", "external_process", counts.get("lanes", 0), pre["trim_parallel_jobs"],
         max(pre["fastqc_threads"], pre["bbduk_threads"]), pre["bbduk_memory_gb"]),
        ("preprocess", "star_and_sort", "external_process", counts.get("lanes", 0), pre["star_parallel_jobs"],
         max(pre["star_threads"], pre["samtools_threads"]),
         max(pre["star_memory_gb"],
             pre["samtools_threads"] * pre["samtools_sort_memory_per_thread_gb"])),
        ("preprocess", "sample_merge", "external_process", counts.get("samples", 0), pre["merge_parallel_jobs"],
         pre["samtools_threads"], pre["merge_memory_gb"]),
        ("dge", "featurecounts", "external_process", 1, 1,
         dge["featurecounts_threads"], dge["featurecounts_memory_gb"]),
        ("dge", "contrast", "external_process", counts.get("contrasts", 0), dge["contrast_parallel_jobs"],
         dge["contrast_threads"], dge["contrast_memory_gb"]),
        ("apa_a", "exact_end_extraction", "python_process", counts.get("samples", 0),
         apa_a["extraction_parallel_jobs"],
         apa_a["extraction_threads"], apa_a["extraction_memory_gb"]),
        ("apa_a", "contrast", "external_process", counts.get("contrasts", 0), apa_a["contrast_parallel_jobs"],
         apa_a["contrast_threads"], apa_a["contrast_memory_gb"]),
        ("apa_b", "engine", "external_process", 1, 1,
         apa_b["engine_threads"], apa_b["engine_memory_gb"]),
        ("apa_b", "contrast", "external_process", counts.get("contrasts", 0), apa_b["contrast_parallel_jobs"],
         apa_b["contrast_threads"], apa_b["contrast_memory_gb"]),
        ("enrichment", "analysis", "external_process", counts.get("enrichment_jobs", counts.get("contrasts", 0)),
         enrichment["parallel_jobs"], enrichment["threads"], enrichment["memory_gb"]),
        ("tracks", "c0_sample", "python_process_and_external_process", counts.get("samples", 0), tracks["parallel_jobs"],
         tracks["samtools_threads"], tracks["memory_gb"]),
        ("tracks", "end_sample", "python_process_and_external_process", counts.get("samples", 0), tracks["parallel_jobs"],
         tracks["samtools_threads"], tracks["memory_gb"]),
    ]
    rows: list[dict[str, Any]] = []
    for stage, unit, executor, units, jobs, threads, memory in definitions:
        effective_jobs = min(jobs, units) if units else jobs
        max_threads = effective_jobs * threads
        max_memory = effective_jobs * memory
        rows.append({
            "stage": stage,
            "work_unit": unit,
            "executor": executor,
            "units": units,
            "max_parallel_jobs": jobs,
            "effective_parallel_jobs": effective_jobs,
            "threads_per_job": threads,
            "memory_gb_per_job": memory,
            "max_threads": max_threads,
            "max_memory_gb": max_memory,
            "total_threads_ceiling": resources["total_threads"],
            "total_memory_gb_ceiling": resources["total_memory_gb"],
            "budget_status": "PASS" if max_threads <= resources["total_threads"] and max_memory <= resources["total_memory_gb"] else "FAIL",
        })
    return rows


def write_resource_plan(resources: dict[str, Any], counts: dict[str, int], path: Path) -> None:
    rows = resource_plan_rows(resources, counts)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def run_bounded(
    stage: str,
    jobs: list[tuple[str, Callable[[], T]]],
    max_workers: int,
    timing_dir: Path,
    progress: Callable[[str, str], None] | None = None,
    on_completed: Callable[[str, T], None] | None = None,
) -> list[T]:
    """Run independent jobs in a bounded pool and return results in input order."""
    if not jobs:
        return []
    timing_dir.mkdir(parents=True, exist_ok=True)
    results: list[T | None] = [None] * len(jobs)
    failures: list[tuple[str, Exception]] = []

    def timed(index: int, label: str, worker: Callable[[], T]) -> tuple[int, T]:
        started = time.time()
        status = "SUCCESS"
        try:
            return index, worker()
        except BaseException:
            status = "FAILED"
            raise
        finally:
            finished = time.time()
            payload = {
                "stage": stage,
                "label": label,
                "started_epoch": started,
                "finished_epoch": finished,
                "elapsed_seconds": round(finished - started, 6),
                "status": status,
            }
            target = timing_dir / f"{label}.json"
            temporary = timing_dir / f".{label}.json.tmp"
            temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
            temporary.replace(target)

    with ThreadPoolExecutor(max_workers=max_workers, thread_name_prefix=stage) as pool:
        futures = {
            pool.submit(timed, index, label, worker): (index, label)
            for index, (label, worker) in enumerate(jobs)
        }
        for future in as_completed(futures):
            _index, label = futures[future]
            try:
                returned_index, value = future.result()
                results[returned_index] = value
                if on_completed is not None:
                    on_completed(label, value)
                if progress is not None:
                    progress(label, "completed")
            except Exception as exc:  # noqa: BLE001 - aggregate independent worker failures
                failures.append((label, exc))
                if progress is not None:
                    progress(label, "failed")
    if failures:
        detail = "; ".join(f"{label}: {exc}" for label, exc in failures)
        raise RuntimeError(f"{stage} worker failures: {detail}")
    return [value for value in results]  # type: ignore[misc]


def _timed_process(
    stage: str,
    index: int,
    label: str,
    worker: Callable[..., T],
    arguments: tuple[Any, ...],
    timing_dir: Path,
) -> tuple[int, T]:
    """Process-pool entry point; all arguments must remain pickleable."""
    started = time.time()
    status = "SUCCESS"
    try:
        return index, worker(*arguments)
    except BaseException:
        status = "FAILED"
        raise
    finally:
        finished = time.time()
        payload = {
            "stage": stage,
            "label": label,
            "executor": "process",
            "started_epoch": started,
            "finished_epoch": finished,
            "elapsed_seconds": round(finished - started, 6),
            "status": status,
            "pid": os.getpid(),
        }
        target = timing_dir / f"{label}.json"
        temporary = timing_dir / f".{label}.{os.getpid()}.json.tmp"
        temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        temporary.replace(target)


def run_bounded_processes(
    stage: str,
    jobs: list[ProcessJob[T]],
    max_workers: int,
    timing_dir: Path,
    progress: Callable[[str, str], None] | None = None,
) -> list[T]:
    """Run CPU-bound independent jobs in processes and preserve input order."""
    if not jobs:
        return []
    timing_dir.mkdir(parents=True, exist_ok=True)
    results: list[T | None] = [None] * len(jobs)
    failures: list[tuple[str, Exception]] = []
    with ProcessPoolExecutor(max_workers=max_workers) as pool:
        futures = {
            pool.submit(_timed_process, stage, index, label, worker, arguments, timing_dir): (index, label)
            for index, (label, worker, arguments) in enumerate(jobs)
        }
        for future in as_completed(futures):
            _index, label = futures[future]
            try:
                returned_index, value = future.result()
                results[returned_index] = value
                if progress is not None:
                    progress(label, "completed")
            except Exception as exc:  # noqa: BLE001 - aggregate independent worker failures
                failures.append((label, exc))
                if progress is not None:
                    progress(label, "failed")
    if failures:
        detail = "; ".join(f"{label}: {exc}" for label, exc in failures)
        raise RuntimeError(f"{stage} worker failures: {detail}")
    return [value for value in results]  # type: ignore[misc]


def merge_tsv_fragments(paths: list[Path], target: Path) -> None:
    """Merge one-row/multi-row TSV fragments deterministically with one header."""
    header: str | None = None
    temporary = target.with_name(f".{target.name}.tmp")
    target.parent.mkdir(parents=True, exist_ok=True)
    with temporary.open("w", encoding="utf-8", newline="") as output:
        for path in paths:
            with path.open(encoding="utf-8", newline="") as source:
                current = source.readline()
                if not current:
                    raise RuntimeError(f"Empty TSV fragment: {path}")
                if header is None:
                    header = current
                    output.write(current)
                elif current != header:
                    raise RuntimeError(f"TSV fragment header mismatch: {path}")
                for line in source:
                    output.write(line)
    temporary.replace(target)
