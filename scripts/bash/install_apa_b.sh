#!/usr/bin/env bash
set -euo pipefail

TAG=""
PREFIX=""
MAMBA="/opt/miniconda/condabin/mamba"
MAMBA_ROOT_PREFIX="/opt/conda_envs/.mamba-rna_ends2tracks-apa-b"

while (($#)); do
  case "$1" in
    --tag) TAG="$2"; shift 2 ;;
    --prefix) PREFIX="$2"; shift 2 ;;
    --mamba) MAMBA="$2"; shift 2 ;;
    -h|--help)
      echo "Usage: install_apa_b.sh --tag RELEASE_TAG [--prefix PATH] [--mamba PATH]"
      exit 0 ;;
    *) echo "Unknown argument: $1" >&2; exit 2 ;;
  esac
done
[[ -n "$TAG" ]] || { echo "--tag is required" >&2; exit 2; }
if [[ -z "$PREFIX" ]]; then
  VERSION="${TAG#v}"
  ENV_TOKEN="${VERSION/-alpha./a}"
  PREFIX="/opt/conda_envs/rna_ends2tracks-apa-b-${ENV_TOKEN}"
fi
[[ ! -e "$PREFIX" ]] || { echo "Target already exists: $PREFIX" >&2; exit 2; }

POLYASEQTRAP_COMMIT="176ea2884ff1c6be7c64bc44fa7661d82d90e718"
DEEPIP_COMMIT="988564875d002b6d5d48d8dfb228cba3492dd776"
HUMAN_MODEL_SHA256="d74138c788102ae57a50664b6858a0b79951430fee9bdcc93b07f9b1ba16edf1"
MOUSE_MODEL_SHA256="aba432c85ef6c14e56a6222106acaffbcc3b9131a86508afdf66311fe57123e9"
SOURCE_DIR="$(mktemp -d "${TMPDIR:-/tmp}/rna_ends2tracks_apa_b.XXXXXX")"
trap 'rm -rf -- "$SOURCE_DIR"' EXIT

git clone --depth 1 --branch "$TAG" https://github.com/MichalGd/rna_ends2tracks.git "$SOURCE_DIR/workflow"
git -C "$SOURCE_DIR/workflow" fetch --depth 1 origin "$TAG"
git -C "$SOURCE_DIR/workflow" checkout --detach "$TAG"

export MAMBA_ROOT_PREFIX
MAMBA_CHANNEL_PRIORITY=strict nice -n 10 "$MAMBA" --no-rc env create \
  --prefix "$PREFIX" --file "$SOURCE_DIR/workflow/environment.apa_b.yml"

"$MAMBA" run -p "$PREFIX" python -c '
import keras
import tensorflow as tf
assert keras.__version__ == "2.10.0", keras.__version__
assert tf.__version__ == "2.10.1", tf.__version__
assert not tf.config.list_physical_devices("GPU")
print("DeepIP TensorFlow CPU runtime: PASS")
'

mkdir -p "$PREFIX/share/rna_ends2tracks-apa-b/sources"
git clone https://github.com/APAexplorer/PolyAseqTrap.git "$SOURCE_DIR/PolyAseqTrap"
git -C "$SOURCE_DIR/PolyAseqTrap" checkout --detach "$POLYASEQTRAP_COMMIT"
git clone https://github.com/APAexplorer/DeepIP.git "$SOURCE_DIR/DeepIP"
git -C "$SOURCE_DIR/DeepIP" checkout --detach "$DEEPIP_COMMIT"

"$MAMBA" run -p "$PREFIX" R CMD INSTALL "$SOURCE_DIR/PolyAseqTrap"
"$MAMBA" run -p "$PREFIX" python -m pip install --no-deps "$SOURCE_DIR/workflow"
cp -a "$SOURCE_DIR/DeepIP/." "$PREFIX/share/rna_ends2tracks-apa-b/sources/DeepIP/"

DEEPIP_ROOT="$PREFIX/share/rna_ends2tracks-apa-b/sources/DeepIP"
HUMAN_MODEL="$(find "$DEEPIP_ROOT/training_model" -type f -iname '*human*.hdf5' -print -quit)"
MOUSE_MODEL="$(find "$DEEPIP_ROOT/training_model" -type f -iname '*mouse*.hdf5' -print -quit)"
[[ -n "$HUMAN_MODEL" && -n "$MOUSE_MODEL" ]] || { echo "DeepIP human or mouse model not found" >&2; exit 1; }
echo "$HUMAN_MODEL_SHA256  $HUMAN_MODEL" | sha256sum --check
echo "$MOUSE_MODEL_SHA256  $MOUSE_MODEL" | sha256sum --check

LOCK="$PREFIX/environment-linux-64.explicit.txt"
"$MAMBA" list -p "$PREFIX" --explicit > "$LOCK"
LOCK_SHA256="$(sha256sum "$LOCK" | awk '{print $1}')"
INSTALLATION_MANIFEST="$PREFIX/installation_manifest.json"
"$PREFIX/bin/python" -c '
import json,sys
path,engine,deepip,script,human,human_sha,mouse,mouse_sha,lock,lock_sha=sys.argv[1:]
payload={
 "schema_version":1,
 "engine":{"name":"PolyAseqTrap","source_commit":engine},
 "deepip":{"name":"DeepIP","source_commit":deepip,"script":script},
 "models":{"human":{"path":human,"sha256":human_sha},"mouse":{"path":mouse,"sha256":mouse_sha}},
 "environment":{"explicit_lock":lock,"sha256":lock_sha},
 "status":"installed_not_pilot_accepted"
}
open(path,"w",encoding="utf-8").write(json.dumps(payload,indent=2,sort_keys=True)+"\n")
' "$INSTALLATION_MANIFEST" "$POLYASEQTRAP_COMMIT" "$DEEPIP_COMMIT" \
  "$DEEPIP_ROOT/DeepIP_test.py" "$HUMAN_MODEL" "$HUMAN_MODEL_SHA256" \
  "$MOUSE_MODEL" "$MOUSE_MODEL_SHA256" "$LOCK" "$LOCK_SHA256"

"$PREFIX/bin/Rscript" -e 'library(PolyAseqTrap); library(Rsamtools); cat("PolyAseqTrap/Rsamtools: PASS\n")'
"$PREFIX/bin/rna-ends2tracks-apa-b" --help >/dev/null
"$PREFIX/bin/rna-ends2tracks-run-apa-b-synthetic-pilot" --help >/dev/null
chmod -R a-w "$PREFIX"
find "$PREFIX" ! -type l -perm /222 -print -quit | grep -q . && {
  echo "Writable content remains in APA-B environment" >&2; exit 1;
}
echo "Installed APA-B candidate environment at $PREFIX"
echo "Pilot acceptance is still required before RUN_APA_B=true."
