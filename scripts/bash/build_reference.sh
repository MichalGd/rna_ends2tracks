#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

usage() {
  echo "Usage: $0 --species human|mouse --assembly GRCh38|GRCm39 --fasta FILE --gtf FILE --output DIR [--threads N]" >&2
}

species=""; assembly=""; fasta=""; gtf=""; output=""; threads=8
while (($#)); do
  case "$1" in
    --species) species="$2"; shift 2 ;;
    --assembly) assembly="$2"; shift 2 ;;
    --fasta) fasta="$2"; shift 2 ;;
    --gtf) gtf="$2"; shift 2 ;;
    --output) output="$2"; shift 2 ;;
    --threads) threads="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) usage; exit 2 ;;
  esac
done
[[ -n "$species" && -n "$assembly" && -n "$fasta" && -n "$gtf" && -n "$output" ]] || { usage; exit 2; }
[[ "$species:$assembly" == "human:GRCh38" || "$species:$assembly" == "mouse:GRCm39" ]] || {
  echo "ERROR: enabled reference profiles are human/GRCh38 and mouse/GRCm39" >&2; exit 2;
}
[[ -f "$fasta" && -f "$gtf" ]] || { echo "ERROR: FASTA or GTF does not exist" >&2; exit 2; }
mkdir -p "$output/STAR_index"
samtools faidx "$fasta"
cut -f1,2 "${fasta}.fai" > "$output/${assembly}.chrom.sizes"
STAR --runMode genomeGenerate --runThreadN "$threads" --genomeDir "$output/STAR_index" \
  --genomeFastaFiles "$fasta" --sjdbGTFfile "$gtf"
sha256sum "$fasta" "$gtf" "$output/${assembly}.chrom.sizes" > "$output/reference_assets.sha256"
echo "Reference built. Copy the matching example manifest, set absolute asset paths and release metadata, then validate it."
