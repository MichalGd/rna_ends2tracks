#!/usr/bin/env bash
set -euo pipefail

usage() {
  echo "Usage: write_versioned_launcher.sh --environment DIR --output FILE"
}

ENV_PREFIX=""
OUTPUT=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --environment) ENV_PREFIX="$2"; shift 2 ;;
    --output) OUTPUT="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; usage >&2; exit 2 ;;
  esac
done

[[ -n "$ENV_PREFIX" ]] || { echo "--environment is required" >&2; exit 2; }
[[ -n "$OUTPUT" ]] || { echo "--output is required" >&2; exit 2; }
ENV_PREFIX="$(cd -- "$ENV_PREFIX" && pwd)"
[[ -x "$ENV_PREFIX/bin/rna-ends2tracks" ]] || {
  echo "Workflow executable is absent: $ENV_PREFIX/bin/rna-ends2tracks" >&2
  exit 2
}
OUTPUT_PARENT="$(dirname -- "$OUTPUT")"
[[ -d "$OUTPUT_PARENT" && -w "$OUTPUT_PARENT" ]] || {
  echo "Launcher directory is absent or not writable: $OUTPUT_PARENT" >&2
  exit 2
}
[[ ! -e "$OUTPUT" ]] || { echo "Refusing to replace launcher: $OUTPUT" >&2; exit 2; }

TEMP="$OUTPUT.tmp.$$"
trap 'rm -f -- "$TEMP"' EXIT
# The dollar expressions below are intentionally written into the generated wrapper.
# shellcheck disable=SC2016
{
  printf '%s\n' '#!/usr/bin/env bash'
  printf '%s\n' 'set -euo pipefail'
  printf 'ENV_PREFIX=%q\n' "$ENV_PREFIX"
  printf '%s\n' 'unset PYTHONHOME PYTHONPATH R_HOME R_LIBS R_LIBS_USER'
  printf '%s\n' 'export PATH="$ENV_PREFIX/bin:${PATH:-/usr/bin:/bin}"'
  printf '%s\n' 'exec "$ENV_PREFIX/bin/rna-ends2tracks" "$@"'
} > "$TEMP"
chmod 0555 "$TEMP"
ln "$TEMP" "$OUTPUT"
rm -f -- "$TEMP"
trap - EXIT
echo "Versioned launcher written: $OUTPUT"
