#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
INSTALLER="$PROJECT_ROOT/scripts/bash/install_release.sh"

# This asserts literal installer source; $ENV_PREFIX must not expand in this test.
# shellcheck disable=SC2016
grep -Fq 'find "$ENV_PREFIX" ! -type l -perm /222 -print -quit' "$INSTALLER"
grep -Fq 'bash scripts/bash/write_versioned_launcher.sh' "$INSTALLER"
# This also asserts literal installer source; $ENV_PREFIX must not expand here.
# shellcheck disable=SC2016
if grep -Fq 'ln -s "$ENV_PREFIX/bin/rna-ends2tracks"' "$INSTALLER"; then
  echo "Installer still creates a non-self-contained launcher symlink" >&2
  exit 1
fi

case "$(uname -s)" in
  MINGW*|MSYS*|CYGWIN*)
    echo "Installer immutability permission test: SKIP (Windows permission emulation)"
    exit 0
    ;;
esac

TEST_ROOT="$(mktemp -d "${TMPDIR:-/tmp}/rna_ends2tracks_immutable_test.XXXXXX")"
trap 'chmod -R u+w "$TEST_ROOT" 2>/dev/null || true; rm -rf -- "$TEST_ROOT"' EXIT

mkdir -p "$TEST_ROOT/environment/bin"
touch "$TEST_ROOT/environment/bin/tool"
ln -s tool "$TEST_ROOT/environment/bin/tool-link"
chmod -R a-w "$TEST_ROOT/environment"

if find "$TEST_ROOT/environment" ! -type l -perm /222 -print -quit | grep -q .; then
  echo "Read-only real content was incorrectly reported writable" >&2
  exit 1
fi

chmod u+w "$TEST_ROOT/environment/bin/tool"
if ! find "$TEST_ROOT/environment" ! -type l -perm /222 -print -quit | grep -q .; then
  echo "Writable real content was not detected" >&2
  exit 1
fi

echo "Installer immutability symlink regression: PASS"
