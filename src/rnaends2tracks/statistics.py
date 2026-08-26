from __future__ import annotations

from pathlib import Path
from typing import Any

from .config import RunPlan, signature_for
from .execution import merge_tsv_fragments, run_bounded
from .external import run
from .receipts import receipt_valid, write_receipt


def r_environment(threads: int, memory_gb: int) -> dict[str, str]:
    value = str(threads)
    return {
        "OMP_NUM_THREADS": value,
        "OPENBLAS_NUM_THREADS": value,
        "MKL_NUM_THREADS": value,
        "VECLIB_MAXIMUM_THREADS": value,
        "NUMEXPR_NUM_THREADS": value,
        "R_MAX_VSIZE": f"{memory_gb}G",
    }


def run_r_contrasts(
    *,
    module: str,
    plan: RunPlan,
    results: Path,
    script: Path,
    common_arguments: list[str],
    outdir: Path,
    log_dir: Path,
    receipt_root: Path,
    index_path: Path,
    parallel_jobs: int,
    threads: int,
    memory_gb: int,
    output_suffixes: list[str],
    signature_inputs: list[str | Path],
    signature_parameters: dict[str, Any],
    dry_run: bool,
    force: bool,
) -> None:
    fragments = outdir / ".index_fragments"
    fragments.mkdir(parents=True, exist_ok=True)
    jobs: list[tuple[str, Any]] = []
    fragment_paths: list[Path] = []
    for contrast in plan.contrasts:
        contrast_id = str(contrast["contrast_id"])
        fragment = fragments / f"{contrast_id}.tsv"
        fragment_paths.append(fragment)

        def worker(contrast=contrast, contrast_id=contrast_id, fragment=fragment):
            receipt_dir = receipt_root / contrast_id
            expected = [outdir / f"{contrast_id}{suffix}" for suffix in output_suffixes]
            expected.append(fragment)
            job_signature = signature_for(signature_inputs, {
                **signature_parameters,
                "module": module,
                "contrast": contrast,
            }) if not dry_run else "dry-run"
            if not force and not dry_run and receipt_valid(receipt_dir, job_signature):
                return fragment
            command = [
                "Rscript", str(script), *common_arguments,
                "--contrast-id", contrast_id, "--index-file", str(fragment),
            ]
            run(
                command, log_dir / f"{contrast_id}.log", dry_run,
                env=r_environment(threads, memory_gb),
            )
            if not dry_run:
                write_receipt(
                    f"{module}_contrast", receipt_dir, job_signature, expected,
                    ["rna-ends2tracks", module, contrast_id],
                )
            return fragment

        jobs.append((contrast_id, worker))
    run_bounded(
        f"{module}_contrasts", jobs, parallel_jobs,
        results / ".checkpoints" / "timings" / module / "contrasts",
    )
    if not dry_run:
        merge_tsv_fragments(fragment_paths, index_path)
