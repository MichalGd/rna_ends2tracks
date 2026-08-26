#!/usr/bin/env bash
set -euo pipefail

usage() {
  echo "Usage: download_pas_atlas_v1_sources.sh --output DIR --accept-ucsc-chain-eula"
}

OUTPUT=""
ACCEPT_UCSC_EULA=false
while [[ $# -gt 0 ]]; do
  case "$1" in
    --output) OUTPUT="$2"; shift 2 ;;
    --accept-ucsc-chain-eula) ACCEPT_UCSC_EULA=true; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; usage >&2; exit 2 ;;
  esac
done

[[ -n "$OUTPUT" ]] || { echo "--output is required" >&2; exit 2; }
$ACCEPT_UCSC_EULA || {
  echo "The mouse conversion requires a UCSC chain file." >&2
  echo "Review https://genome.ucsc.edu/license/EULA.pdf and rerun with --accept-ucsc-chain-eula." >&2
  exit 2
}
command -v curl >/dev/null || { echo "curl is required" >&2; exit 2; }
command -v sha256sum >/dev/null || { echo "sha256sum is required" >&2; exit 2; }

OUTPUT_PARENT="$(dirname "$OUTPUT")"
OUTPUT_NAME="$(basename "$OUTPUT")"
mkdir -p "$OUTPUT_PARENT"
OUTPUT_PARENT="$(cd "$OUTPUT_PARENT" && pwd)"
OUTPUT="$OUTPUT_PARENT/$OUTPUT_NAME"
[[ ! -e "$OUTPUT" ]] || { echo "Refusing to overwrite snapshot directory: $OUTPUT" >&2; exit 2; }
WORK="$OUTPUT.partial.$$"
mkdir "$WORK"
FILES=(
  gencode.v42.polyAs.gtf.gz
  gencode.vM31.polyAs.gtf.gz
  HumanPas.zip
  MousePas.zip
  mm10ToMm39.over.chain.gz
)
PARTIALS=()
COMPLETE=false
cleanup() {
  for partial in "${PARTIALS[@]}"; do
    rm -f -- "$partial"
  done
  if [[ "$COMPLETE" != true && -d "$WORK" && "$WORK" == "$OUTPUT_PARENT/"*.partial.* ]]; then
    rm -rf -- "$WORK"
  fi
}
trap cleanup EXIT

download() {
  local filename="$1"
  local url="$2"
  local partial="$WORK/.${filename}.partial"
  PARTIALS+=("$partial")
  echo "Downloading $filename"
  curl --fail --location --retry 3 --retry-delay 5 --output "$partial" "$url"
  [[ -s "$partial" ]] || { echo "Empty download: $url" >&2; exit 2; }
  mv "$partial" "$WORK/$filename"
}

download gencode.v42.polyAs.gtf.gz \
  https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_42/gencode.v42.polyAs.gtf.gz
download gencode.vM31.polyAs.gtf.gz \
  https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M31/gencode.vM31.polyAs.gtf.gz
download HumanPas.zip \
  https://exon.apps.wistar.org/polya_db/v4/download/4.1/HumanPas.zip
download MousePas.zip \
  https://exon.apps.wistar.org/polya_db/v4/download/4.1/MousePas.zip
download mm10ToMm39.over.chain.gz \
  https://hgdownload.soe.ucsc.edu/goldenPath/mm10/liftOver/mm10ToMm39.over.chain.gz

DOWNLOADED_AT="$(date -u +%FT%TZ)"
{
  printf 'file\turl\trelease\tdownloaded_at_utc\n'
  printf 'gencode.v42.polyAs.gtf.gz\thttps://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_42/gencode.v42.polyAs.gtf.gz\tGENCODE_v42\t%s\n' "$DOWNLOADED_AT"
  printf 'gencode.vM31.polyAs.gtf.gz\thttps://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M31/gencode.vM31.polyAs.gtf.gz\tGENCODE_vM31\t%s\n' "$DOWNLOADED_AT"
  printf 'HumanPas.zip\thttps://exon.apps.wistar.org/polya_db/v4/download/4.1/HumanPas.zip\tPolyA_DB_v4.1\t%s\n' "$DOWNLOADED_AT"
  printf 'MousePas.zip\thttps://exon.apps.wistar.org/polya_db/v4/download/4.1/MousePas.zip\tPolyA_DB_v4.1\t%s\n' "$DOWNLOADED_AT"
  printf 'mm10ToMm39.over.chain.gz\thttps://hgdownload.soe.ucsc.edu/goldenPath/mm10/liftOver/mm10ToMm39.over.chain.gz\tUCSC_mm10_to_mm39\t%s\n' "$DOWNLOADED_AT"
} > "$WORK/download_manifest.tsv"

(
  cd "$WORK"
  sha256sum "${FILES[@]}" download_manifest.tsv > SHA256SUMS.raw
)
chmod -R a-w "$WORK"
mv "$WORK" "$OUTPUT"
COMPLETE=true
echo "PAS atlas v1 source snapshots: PASS"
echo "Snapshot directory: $OUTPUT"
