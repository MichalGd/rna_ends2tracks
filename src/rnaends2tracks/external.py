from __future__ import annotations

import json
import os
import shlex
import shutil
import subprocess
from datetime import datetime, timezone
from pathlib import Path


def require_tools(names: list[str]) -> None:
    missing = [name for name in names if shutil.which(name) is None]
    if missing:
        raise RuntimeError("Required executables are unavailable: " + ", ".join(missing))


def run(command: list[str], log_path: Path, dry_run: bool = False, cwd: Path | None = None) -> None:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    display = shlex.join(command)
    if dry_run:
        with log_path.open("a", encoding="utf-8") as handle:
            handle.write("DRY RUN: " + display + "\n")
        return
    with log_path.open("a", encoding="utf-8") as handle:
        handle.write("COMMAND: " + display + "\n")
        handle.flush()
        subprocess.run(command, cwd=cwd, stdout=handle, stderr=subprocess.STDOUT, check=True)


def event(log_dir: Path, module: str, status: str, message: str) -> None:
    log_dir.mkdir(parents=True, exist_ok=True)
    payload = {
        "time": datetime.now(timezone.utc).isoformat(),
        "module": module,
        "status": status,
        "message": message,
        "pid": os.getpid(),
    }
    with (log_dir / "events.jsonl").open("a", encoding="utf-8") as handle:
        handle.write(json.dumps(payload, sort_keys=True) + "\n")
