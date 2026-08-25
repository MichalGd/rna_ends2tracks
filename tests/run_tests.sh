#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd -- "${SCRIPT_DIR}/.." && pwd)"
export PYTHONPATH="${PROJECT_ROOT}/src${PYTHONPATH:+:${PYTHONPATH}}"
python3 -m unittest discover -s "${SCRIPT_DIR}" -p 'test_*.py' -v
for script in "${PROJECT_ROOT}"/bin/*.sh "${PROJECT_ROOT}"/scripts/bash/*.sh "${PROJECT_ROOT}"/tests/*.sh; do
  bash -n "$script"
done
