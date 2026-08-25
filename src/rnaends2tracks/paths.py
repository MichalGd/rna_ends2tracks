from __future__ import annotations

import sys
from pathlib import Path


def workflow_asset(relative: str) -> Path:
    source = Path(__file__).resolve().parents[2] / relative
    if source.exists():
        return source
    installed = Path(sys.prefix) / "share" / "rna_ends2tracks" / relative
    if installed.exists():
        return installed
    raise RuntimeError(f"Installed workflow asset is missing: {relative}")
