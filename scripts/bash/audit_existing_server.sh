#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

usage() {
  cat >&2 <<'EOF'
Read-only coexistence audit for a server with active cutnrun2tracks jobs.

Usage:
  audit_existing_server.sh --cut-env DIR --output DIR \
    [--human-fasta FILE] [--human-gtf FILE] [--human-chrom-sizes FILE] \
    [--mouse-fasta FILE] [--mouse-gtf FILE] [--mouse-chrom-sizes FILE]

The script does not activate, clone, install, update, index, hash large files,
change permissions, create application directories, or signal running jobs.
EOF
}

cut_env=""
output=""
human_fasta=""; human_gtf=""; human_chrom_sizes=""
mouse_fasta=""; mouse_gtf=""; mouse_chrom_sizes=""
while (($#)); do
  case "$1" in
    --cut-env) cut_env="$2"; shift 2 ;;
    --output) output="$2"; shift 2 ;;
    --human-fasta) human_fasta="$2"; shift 2 ;;
    --human-gtf) human_gtf="$2"; shift 2 ;;
    --human-chrom-sizes) human_chrom_sizes="$2"; shift 2 ;;
    --mouse-fasta) mouse_fasta="$2"; shift 2 ;;
    --mouse-gtf) mouse_gtf="$2"; shift 2 ;;
    --mouse-chrom-sizes) mouse_chrom_sizes="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) usage; exit 2 ;;
  esac
done
[[ -n "$cut_env" && -n "$output" ]] || { usage; exit 2; }
[[ -d "$cut_env" ]] || { echo "ERROR: CUT environment does not exist: $cut_env" >&2; exit 2; }

mkdir -p "$output"
report="$output/coexistence_audit.tsv"
processes="$output/running_cutnrun_processes.txt"
packages="$output/cut_environment_packages.txt"

printf 'section\titem\tstatus\tdetail\n' > "$report"
printf 'system\ttimestamp\tINFO\t%s\n' "$(date --iso-8601=seconds)" >> "$report"
printf 'system\thost\tINFO\t%s\n' "$(hostname -f 2>/dev/null || hostname)" >> "$report"
printf 'system\tload\tINFO\t%s\n' "$(uptime | tr '\t' ' ')" >> "$report"
printf 'system\tdisk\tINFO\t%s\n' "$(df -Pk "$output" | tail -n 1 | tr '\t' ' ')" >> "$report"

ps -eo pid,ppid,user,lstart,etime,stat,args --width 300 | \
  awk 'NR==1 || /cutnrun2tracks|preprocess_batch|align_batch|peakcall_batch|differential_batch|coverage_batch/ {print}' \
  > "$processes"
printf 'jobs\tcutnrun_process_snapshot\tINFO\t%s\n' "$processes" >> "$report"

if command -v conda >/dev/null 2>&1; then
  conda list --prefix "$cut_env" > "$packages"
  printf 'environment\tconda_package_inventory\tOK\t%s\n' "$packages" >> "$report"
else
  printf 'environment\tconda_package_inventory\tWARN\tconda not available in audit shell\n' >> "$report"
fi

commands=(python3 Rscript fastqc multiqc samtools bedtools bamCoverage bedGraphToBigWig STAR bbduk.sh featureCounts)
for command_name in "${commands[@]}"; do
  candidate="$cut_env/bin/$command_name"
  if [[ -x "$candidate" ]]; then
    printf 'tool\t%s\tPRESENT\t%s\n' "$command_name" "$candidate" >> "$report"
  else
    printf 'tool\t%s\tMISSING\tnew rna_ends2tracks environment must provide it\n' "$command_name" >> "$report"
  fi
done

python_imports=(yaml pysam)
for module_name in "${python_imports[@]}"; do
  if "$cut_env/bin/python3" -c "import ${module_name}" >/dev/null 2>&1; then
    printf 'python\t%s\tPRESENT\t%s\n' "$module_name" "$cut_env" >> "$report"
  else
    printf 'python\t%s\tMISSING\tnew environment must provide it\n' "$module_name" >> "$report"
  fi
done

if [[ -x "$cut_env/bin/Rscript" ]]; then
  for package_name in DESeq2 DEXSeq DRIMSeq stageR; do
    if "$cut_env/bin/Rscript" -e "quit(status=ifelse(requireNamespace('${package_name}', quietly=TRUE),0,1))"; then
      printf 'R\t%s\tPRESENT\t%s\n' "$package_name" "$cut_env" >> "$report"
    else
      printf 'R\t%s\tMISSING\tnew environment must provide it\n' "$package_name" >> "$report"
    fi
  done
fi

reference_rows=(
  "human_fasta|$human_fasta"
  "human_gtf|$human_gtf"
  "human_chrom_sizes|$human_chrom_sizes"
  "mouse_fasta|$mouse_fasta"
  "mouse_gtf|$mouse_gtf"
  "mouse_chrom_sizes|$mouse_chrom_sizes"
)
for row in "${reference_rows[@]}"; do
  name="${row%%|*}"; path="${row#*|}"
  [[ -n "$path" ]] || continue
  if [[ -s "$path" && -r "$path" ]]; then
    detail="$(stat -c '%A %U:%G %s_bytes %y %n' "$path")"
    printf 'reference\t%s\tREUSABLE_CANDIDATE\t%s\n' "$name" "$detail" >> "$report"
  else
    printf 'reference\t%s\tFAIL\tmissing, empty, or unreadable: %s\n' "$name" "$path" >> "$report"
  fi
done

cat > "$output/README.txt" <<'EOF'
This directory contains a read-only inventory. PRESENT means only that a file,
executable, import, or package was observed. It does not authorize using or
modifying the cutnrun2tracks environment. FASTA/GTF/chromosome files remain
candidates until assembly, release, contig, checksum, and license checks pass.
EOF

echo "Audit written to $report"
