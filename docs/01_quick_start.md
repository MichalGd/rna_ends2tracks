# Shared-server quick start

This page is for users of the validated shared `biolserv` deployment. The system launcher selects the immutable workflow environment internally. Do not activate Conda and do not export Python, R or tool paths.

## 1. Confirm the launcher

```bash
rna-ends2tracks --version
rna-ends2tracks --help | head -n 25
```

The release documented here is:

```text
0.1.0a11.post2
```

Use the configuration template belonging to the installed release; do not combine a newer executable with an older unreviewed template.

## 2. Create the two project input files

```bash
PROJECT="$HOME/Analysis/my_quantseq_project"
CONFIG_DIR="$PROJECT/config"
RELEASE_ENV="/opt/conda_envs/rna_ends2tracks-0.1.0a11.post2"

mkdir -p "$CONFIG_DIR"
cp "$RELEASE_ENV/share/rna_ends2tracks/config/config.conf" \
  "$CONFIG_DIR/config.conf"
cp "$RELEASE_ENV/share/rna_ends2tracks/config/samplesheet.csv" \
  "$CONFIG_DIR/samplesheet.csv"
```

These are the only per-project input-configuration files. Shared references, FastQ Screen indexes, PAS atlases and the optional APA-B sidecar remain read-only server resources referenced from `config.conf`.

## 3. Prepare `samplesheet.csv`

Keep the header unchanged. One row represents one technical-library/lane FASTQ unit. Rows with the same `sample_id` are merged into one biological sample after lane processing.

- Use absolute FASTQ paths.
- For single-end data, set `library_layout=SE` and leave `fastq_r2` empty.
- For paired-end data, set `library_layout=PE` and provide both FASTQs.
- Use distinct biological replicate identifiers for independent biological samples.
- Use `subject` to identify genuinely paired biological units across conditions.
- Keep all samples and conditions in the same samplesheet.
- Set the correct `genome`, currently `GRCh38` or `GRCm39`.

See the [samplesheet contract](SAMPLESHEET_CONTRACT.md) before encoding technical replicates, paired designs or mixed genomes.

## 4. Edit `config.conf`

At minimum, set:

```text
PROJECT_ID="my_quantseq_project"
SAMPLESHEET="~/Analysis/my_quantseq_project/config/samplesheet.csv"
OUTPUT_DIR="~/Analysis/my_quantseq_project/results"
TMP_DIR="~/Analysis/my_quantseq_project/tmp"
```

For the audited mouse GRCm39 resources on `biolserv`:

```text
MM39_STAR_INDEX="~micgdu/GenomicData/genomicIndices/hsapiens/STAR/STAR-2.7.11b_GRCm39_150bp"
MM39_FASTA="~micgdu/GenomicData/genomesDec2022/GRCm39.primary_assembly.genome.fa"
MM39_GTF="~micgdu/GenomicData/genomesDec2022/gencode.vM31.primary_assembly.annotation.gtf"
MM39_CHROM_SIZES="~micgdu/GenomicData/PAS_atlases/reference_inputs/GRCm39.full.chrom.sizes"
MM39_PAS_ATLAS="~micgdu/GenomicData/PAS_atlases/atlases/GRCm39_gencode_vM31_pas_atlas_v1"
FASTQ_SCREEN_CONFIG="~micgdu/GenomicData/fastq_screen_db/fastq_screen.conf"
FASTQ_SCREEN_MISSING_ACTION="error"
```

Use the audited GRCh38 block for a human project. Unused genome blocks are ignored. Review the global CPU/RAM ceilings and per-stage jobs before starting a large run.

### APA-B on the current server

APA-B is implemented and the `biolserv` GRCm39 QuantSeq REV single-end deployment has passed its synthetic and real-data acceptance pilots. It is off in the portable template so an unvalidated server or protocol cannot enable it accidentally. For an accepted GRCm39 single-end project, use the site-approved installation and validation manifests:

```text
RUN_APA_B=true
APA_B_PILOT_ACCEPTED=true
APA_B_COMMAND_TEMPLATE="auto"
APA_B_INSTALLATION_MANIFEST="/opt/conda_envs/rna_ends2tracks-apa-b-0.1.0a10.post6/installation_manifest.json"
APA_B_VALIDATION_MANIFEST="~micgdu/GenomicData/PAS_atlases/validation/rna_ends2tracks_APA_B_GRCm39_QuantSeq_REV_post6_v2.json"
```

Do not reuse that validation claim for GRCh38 or paired-end APA-B. Those scopes need their own accepted real-data canary and manifest. APA-A remains available independently for supported projects.

## 5. Validate before processing reads

```bash
CONFIG="$CONFIG_DIR/config.conf"

rna-ends2tracks --config "$CONFIG" --stop-after validate
```

Review:

```text
<OUTPUT_DIR>/00_metadata/validated_samples.tsv
<OUTPUT_DIR>/00_metadata/contrasts.tsv
<OUTPUT_DIR>/00_metadata/resource_plan.tsv
<OUTPUT_DIR>/00_metadata/warnings.tsv
```

Confirm sample/replicate counts, paired versus unpaired designs, reference identities, enabled modules and aggregate resource ceilings.

Optionally render every planned command without executing tools:

```bash
rna-ends2tracks --config "$CONFIG" --dry-run
```

## 6. Launch the complete workflow

Use a new `OUTPUT_DIR` for a new scientific run.

```bash
RUN_LOG="$PROJECT/rna_ends2tracks.nohup.log"
PID_FILE="$PROJECT/rna_ends2tracks.pid"

nohup rna-ends2tracks --config "$CONFIG" \
  > "$RUN_LOG" 2>&1 &

RUN_PID=$!
echo "$RUN_PID" > "$PID_FILE"
disown "$RUN_PID"

echo "PID: $RUN_PID"
echo "Log: $RUN_LOG"
```

The external `nohup` log captures launcher output. The chronological workflow log is `<OUTPUT_DIR>/rna_ends2tracks.log`; detailed native-tool logs remain under `<OUTPUT_DIR>/logs/`.

## 7. Monitor progress

```bash
rna-ends2tracks status "$CONFIG"
tail -F "$(sed -n 's/^OUTPUT_DIR=//p' "$CONFIG" | tr -d '"')/rna_ends2tracks.log"
```

For a refreshing status screen:

```bash
watch -n 60 "rna-ends2tracks status '$CONFIG'"
```

The status output reports the active stage, PID state, output counts, free disk space and completed receipts. Bounded job pools also report completed/total units, elapsed time and approximate ETA in the chronological log.

## 8. Verify completion and inspect results

```bash
rna-ends2tracks status "$CONFIG"
```

Require `Workflow status: completed`, then inspect:

```text
<OUTPUT_DIR>/10_reports/report.html
<OUTPUT_DIR>/10_reports/contrast_summary.tsv
<OUTPUT_DIR>/10_reports/differential_gene_expression_summary.tsv
<OUTPUT_DIR>/10_reports/alternative_polyadenylation_summary.tsv
<OUTPUT_DIR>/10_reports/top_enrichment_terms.tsv
<OUTPUT_DIR>/10_reports/ucsc_track_descriptors/
```

Do not interpret differential results before reviewing FastQC, FastQ Screen, alignment/orientation, RSeQC, replicate-level plots and count/distribution diagnostics.

## If the run stops

Read the last status and chronological-log entries first. Correct the underlying cause, then resume from the earliest affected stage, for example:

```bash
nohup rna-ends2tracks \
  --config "$CONFIG" \
  --from-step enrichment \
  > "$PROJECT/rna_ends2tracks.resume-enrichment.log" 2>&1 &
```

Matching receipts reuse validated earlier work. Do not delete or edit receipts to force progress. See [recovery and troubleshooting](recovery_and_troubleshooting.md).
