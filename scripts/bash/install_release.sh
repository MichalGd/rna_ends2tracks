#!/usr/bin/env bash
set -euo pipefail

usage() {
  echo "Usage: install_release.sh --tag TAG [--repo URL] [--env-parent DIR] [--bin-dir DIR] [--mamba EXE] [--no-promote]"
}

TAG=""
REPO="https://github.com/MichalGd/rna_ends2tracks.git"
ENV_PARENT="/opt/conda_envs"
BIN_DIR="/opt/conda_envs/bin"
MAMBA="/opt/miniconda/condabin/mamba"
PROMOTE=true
while [[ $# -gt 0 ]]; do
  case "$1" in
    --tag) TAG="$2"; shift 2 ;;
    --repo) REPO="$2"; shift 2 ;;
    --env-parent) ENV_PARENT="$2"; shift 2 ;;
    --bin-dir) BIN_DIR="$2"; shift 2 ;;
    --mamba) MAMBA="$2"; shift 2 ;;
    --no-promote) PROMOTE=false; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; usage >&2; exit 2 ;;
  esac
done
[[ -n "$TAG" ]] || { echo "--tag is required" >&2; exit 2; }
[[ "$TAG" =~ ^v[0-9]+\.[0-9]+\.[0-9]+-alpha\.[0-9]+([.]post[0-9]+)?$ ]] || { echo "Unexpected release tag: $TAG" >&2; exit 2; }
[[ -x "$MAMBA" ]] || { echo "Mamba is not executable: $MAMBA" >&2; exit 2; }
[[ -d "$ENV_PARENT" && -w "$ENV_PARENT" ]] || { echo "Environment parent is not writable: $ENV_PARENT" >&2; exit 2; }
mkdir -p "$BIN_DIR"
[[ -w "$BIN_DIR" ]] || { echo "Launcher directory is not writable: $BIN_DIR" >&2; exit 2; }

VERSION="${TAG#v}"
ENV_TOKEN="${VERSION/-alpha./a}"
ENV_PREFIX="$ENV_PARENT/rna_ends2tracks-$ENV_TOKEN"
VERSIONED_LAUNCHER="$BIN_DIR/rna-ends2tracks-$VERSION"
STABLE_LAUNCHER="$BIN_DIR/rna-ends2tracks"
[[ ! -e "$ENV_PREFIX" ]] || { echo "Versioned environment already exists: $ENV_PREFIX" >&2; exit 2; }
[[ ! -e "$VERSIONED_LAUNCHER" ]] || { echo "Versioned launcher already exists: $VERSIONED_LAUNCHER" >&2; exit 2; }

STAGING="$(mktemp -d "${TMPDIR:-/tmp}/rna_ends2tracks_install.XXXXXX")"
trap 'rm -rf -- "$STAGING"' EXIT
git clone --quiet --branch "$TAG" --depth 1 "$REPO" "$STAGING/repository"
cd "$STAGING/repository"
COMMIT="$(git rev-parse HEAD)"
TAG_COMMIT="$(git rev-list -n 1 "$TAG")"
[[ "$COMMIT" == "$TAG_COMMIT" ]] || { echo "Tag/checkout commit mismatch" >&2; exit 2; }

INSTALL_AUDIT="$ENV_PARENT/rna_ends2tracks-deployments/$TAG"
mkdir -p "$INSTALL_AUDIT"
cp environment.yml "$INSTALL_AUDIT/environment.yml"
printf 'release=%s\ncommit=%s\nenvironment=%s\n' "$TAG" "$COMMIT" "$ENV_PREFIX" > "$INSTALL_AUDIT/deployment.txt"

export MAMBA_ROOT_PREFIX="$ENV_PARENT/.mamba-rna_ends2tracks"
mkdir -p "$MAMBA_ROOT_PREFIX"
nice -n 10 "$MAMBA" --no-rc env create --prefix "$ENV_PREFIX" --file environment.yml \
  2>&1 | tee "$INSTALL_AUDIT/mamba-create.log"
"$MAMBA" run -p "$ENV_PREFIX" python -m pip install --no-deps . \
  2>&1 | tee "$INSTALL_AUDIT/pip-install.log"
"$MAMBA" run -p "$ENV_PREFIX" python -m unittest discover -s tests -v \
  2>&1 | tee "$INSTALL_AUDIT/python-tests.log"
"$MAMBA" run -p "$ENV_PREFIX" Rscript -e \
  'stopifnot(requireNamespace("DESeq2"),requireNamespace("DEXSeq"),requireNamespace("DRIMSeq"),requireNamespace("stageR"))'
"$MAMBA" run -p "$ENV_PREFIX" Rscript -e \
  'files <- list.files("scripts/R", pattern="[.]R$", full.names=TRUE); stopifnot(length(files) == 4L); invisible(lapply(files, parse))'
"$MAMBA" run -p "$ENV_PREFIX" Rscript scripts/R/dexseq_all_pairs.R --self-test \
  2>&1 | tee "$INSTALL_AUDIT/r-dexseq-hotfix-smoke.log"
"$MAMBA" run -p "$ENV_PREFIX" Rscript tests/R/deseq2_pairing_smoke.R \
  2>&1 | tee "$INSTALL_AUDIT/r-deseq2-smoke.log"
"$ENV_PREFIX/bin/rna-ends2tracks" --version
"$MAMBA" list -p "$ENV_PREFIX" --explicit > "$INSTALL_AUDIT/environment-linux-64.explicit.txt"
"$MAMBA" run -p "$ENV_PREFIX" python -m pip freeze > "$INSTALL_AUDIT/pip-freeze.txt"

bash scripts/bash/write_versioned_launcher.sh \
  --environment "$ENV_PREFIX" \
  --output "$VERSIONED_LAUNCHER"
chmod -R a+rX "$ENV_PREFIX" "$INSTALL_AUDIT"
chmod -R a-w "$ENV_PREFIX"
if find "$ENV_PREFIX" ! -type l -perm /222 -print -quit | grep -q .; then
  echo "Writable content remains in versioned environment" >&2
  exit 2
fi
if $PROMOTE; then
  TEMP_LINK="$BIN_DIR/.rna-ends2tracks.promote.$$"
  ln -s "$VERSIONED_LAUNCHER" "$TEMP_LINK"
  mv -Tf "$TEMP_LINK" "$STABLE_LAUNCHER"
fi
printf 'installed_at=%s\npromoted=%s\n' "$(date --iso-8601=seconds)" "$PROMOTE" >> "$INSTALL_AUDIT/deployment.txt"
echo "Installed $TAG at $ENV_PREFIX"
echo "Versioned launcher: $VERSIONED_LAUNCHER"
$PROMOTE && echo "Stable launcher: $STABLE_LAUNCHER"
