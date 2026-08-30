from __future__ import annotations

import csv
import hashlib
import importlib.metadata
import json
import os
import platform
import shutil
import socket
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from . import __version__


HASH_LIMIT = 1024 * 1024


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _write(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def _receipt_inventory(results: Path) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for path in sorted(results.rglob("run_receipt.json")):
        try:
            value = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, ValueError):
            rows.append({"receipt": path.relative_to(results), "module": "", "status": "INVALID_JSON"})
            continue
        rows.append({
            "receipt": path.relative_to(results).as_posix(), "module": value.get("module", ""),
            "status": "PASS" if value.get("exit_status") == 0 else "FAIL",
            "workflow_version": value.get("workflow_version", ""), "finished_at": value.get("finished_at", ""),
            "host": value.get("host", ""), "platform": value.get("platform", ""),
            "conda_prefix": value.get("environment", {}).get("conda_prefix", ""),
            "command": " ".join(map(str, value.get("command", []))),
            "outputs": len(value.get("outputs", [])), "signature": value.get("signature", ""),
        })
    return rows


def _environment_packages() -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    prefix = Path(os.environ.get("CONDA_PREFIX", ""))
    meta = prefix / "conda-meta"
    if meta.is_dir():
        for path in sorted(meta.glob("*.json")):
            try:
                value = json.loads(path.read_text(encoding="utf-8"))
            except (OSError, ValueError):
                continue
            rows.append({
                "manager": "conda", "package": str(value.get("name", "")),
                "version": str(value.get("version", "")), "build": str(value.get("build", "")),
                "channel": str(value.get("channel", "")),
            })
    if not rows:
        for package in ("rna-ends2tracks", "pysam", "pyyaml"):
            try:
                version = importlib.metadata.version(package)
            except importlib.metadata.PackageNotFoundError:
                continue
            rows.append({"manager": "python", "package": package, "version": version, "build": "", "channel": ""})
    return rows


def _software_versions() -> list[dict[str, str]]:
    commands = {
        "STAR": ["STAR", "--version"], "samtools": ["samtools", "--version"],
        "featureCounts": ["featureCounts", "-v"], "Rscript": ["Rscript", "--version"],
        "bedtools": ["bedtools", "--version"], "multiqc": ["multiqc", "--version"],
        "fastqc": ["fastqc", "--version"], "fastq_screen": ["fastq_screen", "--version"],
        "bowtie2": ["bowtie2", "--version"],
        "RSeQC": ["geneBody_coverage.py", "--version"],
    }
    rows = [{"tool": "rna-ends2tracks", "path": "", "version": __version__, "status": "PASS"},
            {"tool": "python", "path": os.sys.executable, "version": platform.python_version(), "status": "PASS"}]
    for tool, command in commands.items():
        executable = shutil.which(command[0])
        if not executable:
            rows.append({"tool": tool, "path": "", "version": "", "status": "NOT_FOUND"})
            continue
        try:
            completed = subprocess.run(command, capture_output=True, text=True, timeout=15, check=False)
            output = (completed.stdout + "\n" + completed.stderr).strip().splitlines()
            rows.append({"tool": tool, "path": executable, "version": output[0] if output else "",
                         "status": "PASS" if completed.returncode == 0 else f"EXIT_{completed.returncode}"})
        except (OSError, subprocess.TimeoutExpired) as exc:
            rows.append({"tool": tool, "path": executable, "version": str(exc), "status": "UNAVAILABLE"})
    return rows


def _output_manifest(results: Path, provenance_dir: Path) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    excluded = {path.resolve() for path in provenance_dir.glob("*") if path.is_file()}
    for path in sorted(item for item in results.rglob("*") if item.is_file() and not item.is_symlink()):
        if path.resolve() in excluded:
            continue
        stat = path.stat()
        rows.append({
            "path": path.relative_to(results).as_posix(), "size_bytes": stat.st_size,
            "modified_utc": datetime.fromtimestamp(stat.st_mtime, timezone.utc).isoformat(),
            "validation": "sha256" if stat.st_size <= HASH_LIMIT else "size_mtime",
            "sha256": _sha256(path) if stat.st_size <= HASH_LIMIT else "",
        })
    return rows


def generate_provenance_dashboard(plan: Any, results: Path, outdir: Path) -> list[Path]:
    provenance_dir = outdir / "provenance_dashboard"
    provenance_dir.mkdir(parents=True, exist_ok=True)
    summary_rows: list[dict[str, Any]] = [
        {"category": "workflow", "key": "workflow_version", "value": __version__},
        {"category": "workflow", "key": "project_id", "value": plan.project["project_id"]},
        {"category": "runtime", "key": "generated_at", "value": datetime.now(timezone.utc).isoformat()},
        {"category": "runtime", "key": "host", "value": socket.gethostname()},
        {"category": "runtime", "key": "platform", "value": platform.platform()},
        {"category": "runtime", "key": "python", "value": platform.python_version()},
        {"category": "runtime", "key": "conda_prefix", "value": os.environ.get("CONDA_PREFIX", "")},
    ]
    for key in ("_config_path", "_samplesheet_path"):
        path = Path(str(plan.project.get(key, "")))
        summary_rows.append({
            "category": "input", "key": key.removeprefix("_"), "value": str(path),
            "sha256": _sha256(path) if path.is_file() else "MISSING",
        })
    for assembly, reference in sorted(plan.references.items()):
        for key in (
            "species", "release", "fasta", "gtf", "star_index", "chrom_sizes",
            "pas_atlas", "rseqc_bed",
        ):
            summary_rows.append({"category": f"reference:{assembly}", "key": key, "value": reference.get(key, "")})
        for name, digest in sorted(reference.get("pas_atlas_checksums", {}).items()):
            summary_rows.append({"category": f"reference:{assembly}", "key": f"pas_atlas_sha256:{name}", "value": digest})
    apa_b = plan.project.get("apa_b", {})
    summary_rows.extend([
        {"category": "apa_b", "key": "enabled", "value": apa_b.get("enabled", False)},
        {"category": "apa_b", "key": "interpretation_status",
         "value": "VALIDATED_PILOT_ACCEPTED" if apa_b.get("enabled", False) else "DISABLED_NOT_VALIDATED"},
        {"category": "apa_b", "key": "validation_manifest", "value": apa_b.get("validation_manifest", "")},
    ])
    summary = provenance_dir / "provenance_summary.tsv"
    _write(summary, summary_rows, ["category", "key", "value", "sha256"])

    receipts = _receipt_inventory(results)
    receipt_path = provenance_dir / "receipt_inventory.tsv"
    _write(receipt_path, receipts, [
        "receipt", "module", "status", "workflow_version", "finished_at", "host", "platform",
        "conda_prefix", "command", "outputs", "signature",
    ])
    packages = _environment_packages()
    package_path = provenance_dir / "environment_packages.tsv"
    _write(package_path, packages, ["manager", "package", "version", "build", "channel"])
    software = _software_versions()
    software_path = provenance_dir / "software_versions.tsv"
    _write(software_path, software, ["tool", "path", "version", "status"])
    outputs = _output_manifest(results, provenance_dir)
    output_path = provenance_dir / "output_manifest.tsv"
    _write(output_path, outputs, ["path", "size_bytes", "modified_utc", "validation", "sha256"])
    dashboard = provenance_dir / "dashboard.json"
    dashboard.write_text(json.dumps({
        "schema_version": 1, "workflow_version": __version__, "project_id": plan.project["project_id"],
        "receipt_counts": {
            "total": len(receipts), "pass": sum(row.get("status") == "PASS" for row in receipts),
            "fail": sum(row.get("status") == "FAIL" for row in receipts),
        },
        "output_files": len(outputs), "environment_packages": len(packages), "software_tools": len(software),
        "apa_b_interpretation_status": "VALIDATED_PILOT_ACCEPTED" if apa_b.get("enabled", False) else "DISABLED_NOT_VALIDATED",
    }, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return [summary, receipt_path, package_path, software_path, output_path, dashboard]
