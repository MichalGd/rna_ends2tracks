from __future__ import annotations

import json
import os
import socket
from collections.abc import Iterator
from contextlib import contextmanager
from datetime import datetime, timezone
from pathlib import Path


@contextmanager
def run_lock(results: Path) -> Iterator[Path]:
    """Prevent two workflow processes from mutating one output directory."""
    checkpoint = results / ".checkpoints"
    checkpoint.mkdir(parents=True, exist_ok=True)
    path = checkpoint / "workflow.lock"
    payload = {
        "pid": os.getpid(), "host": socket.gethostname(),
        "started_at": datetime.now(timezone.utc).isoformat(),
        "results": str(results.resolve()),
    }
    try:
        descriptor = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o644)
    except FileExistsError as exc:
        try:
            owner = path.read_text(encoding="utf-8").strip()
        except OSError:
            owner = "unreadable lock metadata"
        raise RuntimeError(
            f"Output directory is locked by another workflow process: {path}\n{owner}\n"
            "Verify that process has ended before removing a stale lock."
        ) from exc
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
            json.dump(payload, handle, indent=2, sort_keys=True)
            handle.write("\n")
        yield path
    finally:
        try:
            path.unlink()
        except FileNotFoundError:
            pass
