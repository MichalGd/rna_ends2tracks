from __future__ import annotations

import csv
import json
import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from .config import RunPlan, signature_for, workflow_requirements
from .external import event
from .receipts import receipt_valid, sha256, write_receipt


SUPPORTED_CLEANUP_EVIDENCE = re.compile(r"^0[.]1[.]0a(?:9|10)(?:[.]post[0-9]+)?$")


def _safe_results_root(results: Path) -> Path:
    root = results.expanduser().resolve()
    anchor = Path(root.anchor).resolve()
    home = Path.home().resolve()
    if root == anchor or root == home or len(root.parts) < 3:
        raise RuntimeError(f"Refusing cleanup for unsafe results directory: {root}")
    return root


def _prior_cleanup_removed_paths(root: Path) -> set[Path]:
    """Return paths proven removed by a previously successful cleanup.

    A resumed run may validate an upstream receipt after its dispensable outputs
    were already removed by an earlier successful cleanup.  Trust only the
    cumulative manifest covered by a successful cleanup receipt; an arbitrary
    manifest is not sufficient evidence.
    """
    cleanup_dir = root / "provenance" / "cleanup"
    receipt_path = cleanup_dir / "run_receipt.json"
    manifest = cleanup_dir / "cleanup_manifest.tsv"
    try:
        receipt: dict[str, Any] = json.loads(receipt_path.read_text(encoding="utf-8"))
        output = next(
            item for item in receipt.get("outputs", [])
            if Path(str(item.get("path", ""))).resolve() == manifest.resolve()
        )
        stat = manifest.stat()
        if (
            receipt.get("exit_status") != 0 or receipt.get("schema_version") != 1
            or receipt.get("module") != "cleanup"
        ):
            return set()
        if stat.st_size != int(output["size"]):
            return set()
        if output.get("validation") == "size_mtime":
            if stat.st_mtime_ns != int(output["mtime_ns"]):
                return set()
        elif sha256(manifest) != str(output["sha256"]):
            return set()
    except (OSError, ValueError, KeyError, TypeError, AttributeError, StopIteration):
        return set()

    removed: set[Path] = set()
    try:
        with manifest.open(encoding="utf-8", newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                if row.get("status") != "removed":
                    continue
                path = Path(row.get("path", "")).resolve(strict=False)
                if path.is_relative_to(root):
                    removed.add(path)
    except (OSError, ValueError):
        return set()
    return removed


def _successful_receipt(module_dir: Path, removed_paths: set[Path]) -> bool:
    receipt_path = module_dir / "run_receipt.json"
    try:
        receipt: dict[str, Any] = json.loads(receipt_path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return False
    # Cleanup does not reuse a computation.  Its guard only needs auditable
    # evidence that the stage completed successfully.  Requiring an exact
    # patch-version match made safe cross-version resumes impossible.
    if (
        receipt.get("schema_version") != 1 or receipt.get("exit_status") != 0
        or not SUPPORTED_CLEANUP_EVIDENCE.fullmatch(
            str(receipt.get("workflow_version", "")).strip()
        )
        or not str(receipt.get("module", "")).strip()
    ):
        return False
    outputs = receipt.get("outputs")
    if not isinstance(outputs, list) or not outputs:
        return False
    for output in outputs:
        path = Path(str(output.get("path", "")))
        try:
            if not path.is_file():
                if path.resolve(strict=False) in removed_paths:
                    continue
                return False
            stat = path.stat()
            if stat.st_size != int(output["size"]):
                return False
            if output.get("validation") == "size_mtime":
                if stat.st_mtime_ns != int(output["mtime_ns"]):
                    return False
            elif sha256(path) != str(output["sha256"]):
                return False
        except (OSError, KeyError, TypeError, ValueError, AttributeError):
            return False
    return True


def _require_successful_workflow(plan: RunPlan, root: Path) -> None:
    removed_paths = _prior_cleanup_removed_paths(root)
    required = ["02_alignment", "10_reports"]
    modules = plan.project.get("modules", {})
    requirements = workflow_requirements(plan.project)
    if modules.get("rseqc", False):
        required.append("01_qc/rseqc")
    if requirements["exact_ends"]:
        required.append("03_exact_ends")
    if requirements["active_pas"]:
        required.append("04_active_pas")
    if modules.get("gene_expression", True):
        required.append("05_gene_expression")
    if modules.get("apa_a", True):
        required.append("06_apa_a_mcell2019")
    if plan.project.get("apa_b", {}).get("enabled", False):
        required.append("07_apa_b")
    if requirements["apa_comparison"]:
        required.append("08_apa_comparison")
    if modules.get("tracks", True):
        required.append("09_tracks")
    incomplete = [
        name for name in required
        if not _successful_receipt(root / name, removed_paths)
    ]
    enrichment_enabled = bool(
        (modules.get("dge_enrichment", False) and modules.get("gene_expression", True))
        or (modules.get("apa_enrichment", False) and (
            modules.get("apa_a", True) or plan.project.get("apa_b", {}).get("enabled", False)
        ))
    )
    if enrichment_enabled and not _successful_receipt(
        root / "10_reports" / "enrichment_summary", removed_paths
    ):
        incomplete.append("10_reports/enrichment_summary")
    if incomplete:
        raise RuntimeError(
            "Cleanup requires a successful complete workflow; missing or invalid receipts: "
            + ", ".join(incomplete)
        )


def _add_files(
    targets: dict[Path, str], root: Path, category: str, directory: Path, patterns: tuple[str, ...]
) -> None:
    if not directory.is_dir():
        return
    for pattern in patterns:
        for path in directory.rglob(pattern):
            if path.is_file() or path.is_symlink():
                resolved = path.resolve()
                if not resolved.is_relative_to(root):
                    raise RuntimeError(f"Cleanup target escapes results directory: {path} -> {resolved}")
                targets[path] = category


def cleanup_targets(plan: RunPlan, results: Path) -> list[tuple[str, Path]]:
    root = _safe_results_root(results)
    settings = plan.project["cleanup"]
    targets: dict[Path, str] = {}

    if not settings["keep_trimmed_fastq"]:
        _add_files(
            targets, root, "trimmed_fastq", root / "01_qc" / "trimmed_fastq",
            ("*.trimmed.fastq.gz", ".*.tmp.fastq.gz"),
        )

    if not settings["keep_lane_bams"]:
        for lane in plan.sample_rows:
            token = f"{lane['sample_id']}.{lane['technical_replicate_id']}.{lane['lane_id']}"
            lane_dir = root / "02_alignment" / lane["sample_id"] / "lanes"
            _add_files(
                targets,
                root,
                "lane_alignment_bam",
                lane_dir,
                (
                    f"{token}.bam",
                    f"{token}.bam.bai",
                    f"{token}.star.Aligned.out.bam",
                    f"{token}.star.Aligned.out.bam.bai",
                    f".{token}.bam.tmp",
                ),
            )
        for sample in plan.samples:
            _add_files(
                targets, root, "merged_all_alignment_bam", root / "02_alignment" / sample["sample_id"],
                (f"{sample['sample_id']}.all_alignments.bam", f"{sample['sample_id']}.all_alignments.bam.bai"),
            )

    track_dir = root / "09_tracks" / ".intermediate"
    if not settings["keep_track_strand_bams"]:
        _add_files(targets, root, "track_strand_bam", track_dir, ("*.strand.bam",))
    if not settings["keep_track_bedgraphs"]:
        _add_files(targets, root, "track_bedgraph", track_dir, ("*.bedGraph",))

    return sorted(((category, path) for path, category in targets.items()), key=lambda item: str(item[1]))


def clean_intermediates(
    plan: RunPlan,
    results: Path,
    dry_run: bool = False,
    force: bool = False,
    require_success: bool = True,
) -> Path | None:
    root = _safe_results_root(results)
    settings = plan.project["cleanup"]
    log_dir = root / "logs"
    if not settings["enabled"]:
        event(log_dir, "cleanup", "disabled", "Intermediate cleanup disabled in project configuration")
        return None
    if require_success:
        _require_successful_workflow(plan, root)

    cleanup_dir = root / "provenance" / "cleanup"
    manifest = cleanup_dir / "cleanup_manifest.tsv"
    signature = signature_for([], {
        "module": "cleanup",
        "settings": settings,
        "samples": [sample["sample_id"] for sample in plan.samples],
        "results": str(root),
    })
    targets = cleanup_targets(plan, root)
    if not force and not targets and receipt_valid(cleanup_dir, signature):
        event(log_dir, "cleanup", "skipped", "No configured intermediates remain")
        return manifest
    if dry_run:
        event(log_dir, "cleanup", "dry_run", f"Would remove {len(targets)} intermediate files")
        return None

    cleanup_dir.mkdir(parents=True, exist_ok=True)
    removed_at = datetime.now(timezone.utc).isoformat()
    rows: list[dict[str, str | int]] = []
    if manifest.is_file():
        with manifest.open(encoding="utf-8", newline="") as handle:
            rows.extend(csv.DictReader(handle, delimiter="\t"))
    removed_bytes = 0
    for category, path in targets:
        size = path.lstat().st_size
        path.unlink()
        removed_bytes += size
        rows.append({
            "category": category,
            "path": str(path.resolve(strict=False)),
            "size_bytes": size,
            "status": "removed",
            "removed_at": removed_at,
        })

    # Remove only now-empty workflow-owned intermediate directories. Final
    # deliverables, receipts and audit records are never targeted.
    empty_directories = [
        root / "01_qc" / "trimmed_fastq",
        root / "09_tracks" / ".intermediate",
    ]
    empty_directories.extend(
        root / "09_tracks" / ".intermediate" / sample["sample_id"]
        for sample in plan.samples
    )
    empty_directories.extend(
        root / "02_alignment" / sample["sample_id"] / "lanes"
        for sample in plan.samples
    )
    for directory in sorted(empty_directories, key=lambda item: len(item.parts), reverse=True):
        if directory.is_dir() and not any(directory.iterdir()):
            directory.rmdir()

    temporary = cleanup_dir / ".cleanup_manifest.tsv.tmp"
    with temporary.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["category", "path", "size_bytes", "status", "removed_at"],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(manifest)
    write_receipt("cleanup", cleanup_dir, signature, [manifest], ["rna-ends2tracks", "cleanup"])
    event(
        log_dir,
        "cleanup",
        "completed",
        f"Removed {len(targets)} dispensable intermediate files "
        f"({removed_bytes} bytes)",
    )
    return manifest
