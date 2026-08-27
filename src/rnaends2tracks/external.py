from __future__ import annotations

import json
import os
import shlex
import shutil
import subprocess
import threading
from datetime import datetime, timezone
from pathlib import Path


_EVENT_LOCK = threading.Lock()


def _results_root(log_dir: Path) -> Path:
    """Return the workflow result root for a logs directory or its descendant."""
    for candidate in (log_dir, *log_dir.parents):
        if candidate.name == "logs":
            return candidate.parent
    return log_dir.parent


def _atomic_json(path: Path, payload: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    temporary.replace(path)


def _update_status(results: Path, payload: dict[str, object]) -> None:
    path = results / "00_metadata" / "run_status.json"
    try:
        state = json.loads(path.read_text(encoding="utf-8")) if path.is_file() else {}
    except (OSError, ValueError):
        state = {}
    stages = state.setdefault("stages", {})
    if not isinstance(stages, dict):
        stages = {}
        state["stages"] = stages
    module = str(payload["module"])
    status = str(payload["status"])
    stages[module] = {
        "status": status,
        "time": payload["time"],
        "message": payload["message"],
        "pid": payload["pid"],
    }
    state.update({
        "current_stage": module,
        "last_status": status,
        "last_message": payload["message"],
        "updated_at": payload["time"],
        "pid": payload["pid"],
    })
    if module == "workflow":
        state["workflow_status"] = status
    elif status == "failed":
        state["workflow_status"] = "failed"
    elif state.get("workflow_status") not in {"completed", "failed"}:
        state["workflow_status"] = "running"
    _atomic_json(path, state)


def _command_context(log_path: Path) -> tuple[Path | None, str, str]:
    for candidate in (log_path.parent, *log_path.parents):
        if candidate.name == "logs":
            relative = log_path.relative_to(candidate)
            module = relative.parts[0] if len(relative.parts) > 1 else "external"
            return candidate, module, log_path.stem
    return None, "external", log_path.stem


def require_tools(names: list[str]) -> None:
    missing = [name for name in names if shutil.which(name) is None]
    if missing:
        raise RuntimeError("Required executables are unavailable: " + ", ".join(missing))


def run(
    command: list[str], log_path: Path, dry_run: bool = False, cwd: Path | None = None,
    env: dict[str, str] | None = None,
) -> None:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    display = shlex.join(command)
    workflow_logs, module, label = _command_context(log_path)
    if dry_run:
        with log_path.open("a", encoding="utf-8") as handle:
            handle.write("DRY RUN: " + display + "\n")
        if workflow_logs is not None:
            event(workflow_logs, module, "dry_run", f"{label}: {display}")
        return
    if workflow_logs is not None:
        event(workflow_logs, module, "started", f"{label}: {display}; details={log_path}")
    try:
        with log_path.open("a", encoding="utf-8") as handle:
            handle.write("COMMAND: " + display + "\n")
            handle.flush()
            environment = None if env is None else {**os.environ, **env}
            subprocess.run(command, cwd=cwd, env=environment, stdout=handle, stderr=subprocess.STDOUT, check=True)
    except subprocess.CalledProcessError as exc:
        if workflow_logs is not None:
            event(workflow_logs, module, "failed",
                  f"{label}: exit_status={exc.returncode}; details={log_path}")
        raise
    if workflow_logs is not None:
        event(workflow_logs, module, "completed", f"{label}: exit_status=0; details={log_path}")


def event(log_dir: Path, module: str, status: str, message: str) -> None:
    log_dir.mkdir(parents=True, exist_ok=True)
    payload = {
        "time": datetime.now(timezone.utc).isoformat(),
        "module": module,
        "status": status,
        "message": message,
        "pid": os.getpid(),
    }
    results = _results_root(log_dir)
    severity = "ERROR" if status == "failed" else "INFO"
    readable = (
        f"{payload['time']} {severity} [{module}] {status.upper()} "
        f"{message} (pid={payload['pid']})\n"
    )
    with _EVENT_LOCK:
        with (log_dir / "events.jsonl").open("a", encoding="utf-8") as handle:
            handle.write(json.dumps(payload, sort_keys=True) + "\n")
        with (results / "rna_ends2tracks.log").open("a", encoding="utf-8") as handle:
            handle.write(readable)
        _update_status(results, payload)


def read_run_status(results: Path) -> dict[str, object]:
    path = results / "00_metadata" / "run_status.json"
    if not path.is_file():
        raise RuntimeError(f"Run status is unavailable: {path}")
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError) as exc:
        raise RuntimeError(f"Cannot read run status: {path}") from exc
    if not isinstance(payload, dict):
        raise RuntimeError(f"Invalid run status document: {path}")
    return payload
