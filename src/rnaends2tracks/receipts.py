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
    receipt = {
        "schema_version": 1,
        "module": module,
        "workflow_version": "0.1.0",
        "signature": signature,
        "finished_at": datetime.now(timezone.utc).isoformat(),
        "host": socket.gethostname(),
        "platform": platform.platform(),
        "environment": {"conda_prefix": os.environ.get("CONDA_PREFIX", ""), "container": os.environ.get("APPTAINER_CONTAINER", "")},
        "software": software,
        "command": command,
        "exit_status": 0,
        "outputs": [
            {"path": str(path.resolve()), "size": path.stat().st_size, "sha256": sha256(path)}
            for path in outputs if path.is_file()
        ],
    }
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
    if receipt.get("signature") != signature or receipt.get("exit_status") != 0:
        return False
    for output in receipt.get("outputs", []):
        path = Path(output["path"])
        if not path.is_file() or path.stat().st_size != output["size"] or sha256(path) != output["sha256"]:
            return False
    return True
