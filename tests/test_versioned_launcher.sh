#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
WRITER="$PROJECT_ROOT/scripts/bash/write_versioned_launcher.sh"
TEST_ROOT="$(mktemp -d "${TMPDIR:-/tmp}/rna_ends2tracks_launcher_test.XXXXXX")"
trap 'chmod -R u+w "$TEST_ROOT" 2>/dev/null || true; rm -rf -- "$TEST_ROOT"' EXIT

ENV_PREFIX="$TEST_ROOT/environment"
LAUNCHER="$TEST_ROOT/bin/rna-ends2tracks-test"
mkdir -p "$ENV_PREFIX/bin" "$TEST_ROOT/bin"
# The dollar expressions below belong to the fake executable, not this test process.
# shellcheck disable=SC2016
{
  printf '%s\n' '#!/usr/bin/env bash'
  printf '%s\n' 'printf "PATH=%s\n" "$PATH"'
  printf '%s\n' 'printf "R_HOME=%s\n" "${R_HOME-unset}"'
  printf '%s\n' 'printf "PYTHONPATH=%s\n" "${PYTHONPATH-unset}"'
  printf '%s\n' 'printf "STAR=%s\n" "$(command -v STAR)"'
  printf '%s\n' 'printf "ARGC=%s\n" "$#"'
  printf '%s\n' 'printf "ARG1=%s\n" "$1"'
} > "$ENV_PREFIX/bin/rna-ends2tracks"
printf '%s\n' '#!/usr/bin/env bash' 'exit 0' > "$ENV_PREFIX/bin/STAR"
chmod 0755 "$ENV_PREFIX/bin/rna-ends2tracks" "$ENV_PREFIX/bin/STAR"

bash "$WRITER" --environment "$ENV_PREFIX" --output "$LAUNCHER"
[[ -f "$LAUNCHER" && ! -L "$LAUNCHER" ]]
if find "$LAUNCHER" ! -type l -perm /222 -print -quit | grep -q .; then
  echo "Versioned launcher must be immutable" >&2
  exit 1
fi

RESULT="$(env -i PATH=/usr/bin:/bin R_HOME=/contaminated PYTHONPATH=/contaminated \
  "$LAUNCHER" "alpha beta")"
grep -Fqx "PATH=$ENV_PREFIX/bin:/usr/bin:/bin" <<< "$RESULT"
grep -Fqx 'R_HOME=unset' <<< "$RESULT"
grep -Fqx 'PYTHONPATH=unset' <<< "$RESULT"
grep -Fqx "STAR=$ENV_PREFIX/bin/STAR" <<< "$RESULT"
grep -Fqx 'ARGC=1' <<< "$RESULT"
grep -Fqx 'ARG1=alpha beta' <<< "$RESULT"

echo "Self-contained versioned launcher regression: PASS"
