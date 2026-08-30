from __future__ import annotations

import hashlib
import json
import os
import platform
import shutil
import socket
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from . import __version__

HASH_LIMIT_BYTES = 64 * 1024 * 1024
COMPATIBLE_WORKFLOW_VERSIONS = {
    "0.1.0a9.post1": frozenset({"0.1.0a9"}),
    "0.1.0a9.post2": frozenset({"0.1.0a9", "0.1.0a9.post1"}),
    "0.1.0a11.post1": frozenset({"0.1.0a11"}),
    "0.1.0a11.post2": frozenset({"0.1.0a11", "0.1.0a11.post1"}),
}


def workflow_version_compatible(receipt_version: object) -> bool:
    """Allow only explicitly audited base receipts for this post-release hotfix."""
    value = str(receipt_version)
    return value == __version__ or value in COMPATIBLE_WORKFLOW_VERSIONS.get(__version__, ())


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_receipt(module: str, output_dir: Path, signature: str, outputs: list[Path], command: list[str]) -> Path:
    missing = [str(path) for path in outputs if not path.exists()]
    if missing:
        raise RuntimeError("Cannot publish receipt; outputs are missing: " + ", ".join(missing))
    software = {}
    for executable in ("python3", "STAR", "samtools", "featureCounts", "Rscript", "bedtools"):
        location = shutil.which(executable)
        if location:
            software[executable] = {"path": location}
    output_records = []
    for path in outputs:
        if not path.is_file():
            continue
        stat = path.stat()
        record: dict[str, Any] = {
            "path": str(path.resolve()), "size": stat.st_size, "mtime_ns": stat.st_mtime_ns,
            "validation": "sha256" if stat.st_size <= HASH_LIMIT_BYTES else "size_mtime",
        }
        if stat.st_size <= HASH_LIMIT_BYTES:
            record["sha256"] = sha256(path)
        output_records.append(record)
    receipt = {
        "schema_version": 1,
        "module": module,
        "workflow_version": __version__,
        "signature": signature,
        "finished_at": datetime.now(timezone.utc).isoformat(),
        "host": socket.gethostname(),
        "platform": platform.platform(),
        "environment": {"conda_prefix": os.environ.get("CONDA_PREFIX", ""), "container": os.environ.get("APPTAINER_CONTAINER", "")},
        "software": software,
        "command": command,
        "exit_status": 0,
        "outputs": output_records,
    }
    output_dir.mkdir(parents=True, exist_ok=True)
    receipt_path = output_dir / "run_receipt.json"
    temporary = output_dir / ".run_receipt.json.tmp"
    temporary.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    temporary.replace(receipt_path)
    return receipt_path


def receipt_valid(output_dir: Path, signature: str) -> bool:
    receipt_path = output_dir / "run_receipt.json"
    try:
        receipt: dict[str, Any] = json.loads(receipt_path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return False
    if (
        not workflow_version_compatible(receipt.get("workflow_version"))
        or receipt.get("signature") != signature
        or receipt.get("exit_status") != 0
    ):
        return False
    for output in receipt.get("outputs", []):
        path = Path(output["path"])
        if not path.is_file():
            return False
        stat = path.stat()
        if stat.st_size != output["size"]:
            return False
        if output.get("validation") == "size_mtime":
            if stat.st_mtime_ns != output.get("mtime_ns"):
                return False
        elif sha256(path) != output.get("sha256"):
            return False
    return True
