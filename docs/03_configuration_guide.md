# Configuration guide

`config.conf` and `samplesheet.csv` are the only per-project configuration
files. The samplesheet describes the biological and sequencing units;
`config.conf` selects the workflow behavior, references and resource limits.
Start with the installed templates rather than copying settings from an older
analysis. See the [quick start](01_quick_start.md) for the copy commands.

## Syntax and paths

The file uses a restricted `KEY=value` grammar. Values may be quoted, but shell
commands, substitutions and arbitrary Bash are rejected. Use absolute paths for
the samplesheet, output directory, temporary directory and reference assets.

```text
PROJECT_ID="my_quantseq_project"
SAMPLESHEET="~/Analysis/my_project/config/samplesheet.csv"
OUTPUT_DIR="~/Analysis/my_project/results"
TMP_DIR="~/Analysis/my_project/tmp"
```

Use a new `OUTPUT_DIR` when samples are added or removed. Active-PAS discovery
is project-wide and condition-blind, so changing the sample set changes the PAS
universe.

## Library layout and protocol

Choose one supported project-wide layout and matching protocol:

| Data | `LIBRARY_LAYOUT` | `LIBRARY_PROTOCOL` | Samplesheet FASTQs |
|---|---|---|---|
| QuantSeq REV V1 single-end | `single_end` | `quantseq_rev_v1_se` | R1 only |
| QuantSeq REV V2 single-end | `single_end` | `quantseq_rev_v2_se` | R1 only |
| QuantSeq REV V1 paired-end | `paired_end` | `quantseq_rev_v1_pe` | R1 and R2 |
| QuantSeq REV V2 paired-end | `paired_end` | `quantseq_rev_v2_pe` | R1 and R2 |

For QuantSeq REV paired-end data, both mates contribute to alignment and
conventional coverage, while R1 is the end-defining mate for APA. The default
`PE_R2_TRIM_5P=12` removes the V2 random-primer-derived R2 prefix. Do not change
it without protocol-specific evidence. The current workflow assumes no UMI:

```text
UMI_PRESENT=false
MAPPING_POLICY="unique_primary"
END_SOFT_CLIP_POLICY="exclude_and_report"
```

Replicate and pairing metadata belong in the samplesheet. `PAIRING_MODE="auto"`
uses `PAIRING_COLUMN="subject"` only for contrasts with complete biological
matching and resolves other contrasts as unpaired. See the [samplesheet
contract](SAMPLESHEET_CONTRACT.md).

## Modules

The main analysis switches are independent:

```text
RUN_GENE_EXPRESSION=true
RUN_APA_A_MCELL2019=true
RUN_APA_A2=true
RUN_APA_B=true
APA_B_PILOT_ACCEPTED=true
RUN_DGE_ENRICHMENT=true
RUN_APA_ENRICHMENT=true
RUN_TRACKS=true
RUN_RSEQC=true
RUN_FASTQ_SCREEN=true
```

APA-A is the preserved Mcell2019-style method. APA-A2 independently reruns
DEXSeq on the same condition-blind C3 catalog and corrects PAU effects, paired
averaging, effect qualification, and transcript-coordinate shifts. APA-B uses
the independently pinned PolyAseqTrap/DeepIP sidecar and can run only for an assembly and protocol
covered by accepted manifests. The audited shared-server GRCh38/GRCm39
QuantSeq REV V2 single-end scope is accepted, so all three APA methods are on in
the new-project template; consult the [quick start](01_quick_start.md#apa-b-on-the-current-server).
Set `RUN_APA_A_MCELL2019=false`, `RUN_APA_A2=false`, or `RUN_APA_B=false` to
disable each method independently.
An enabled APA-B run outside the accepted scope fails validation instead of
silently skipping or borrowing APA-A results.

Enrichment is enabled by default for DGE and APA and can use ORA, GSEA, GO,
Reactome, Hallmark and KEGG collections. Plot and network limits are adjustable
without changing the underlying differential analysis.

## Assembly-matched references

Each samplesheet `genome` selects one complete reference block. The current
validated profiles are GRCh38/GENCODE v42 and GRCm39/GENCODE vM31. Keep the
FASTA, GTF, STAR index, chromosome sizes, PAS atlas and optional RSeQC BED from
the same assembly/release profile.

```text
MM39_STAR_INDEX="/absolute/path/to/GRCm39/STAR_index"
MM39_FASTA="/absolute/path/to/GRCm39.genome.fa"
MM39_GTF="/absolute/path/to/gencode.vM31.gtf"
MM39_CHROM_SIZES="/absolute/path/to/GRCm39.chrom.sizes"
MM39_PAS_ATLAS="/absolute/path/to/GRCm39_gencode_vM31_pas_atlas_v1"
MM39_RSEQC_BED=""
```

An empty RSeQC BED makes the workflow generate an assembly-matched BED12 file
from the configured GTF. Unused genome blocks are ignored. Cross-genome
contrasts are rejected.

FastQ Screen uses its own site configuration containing the indexed comparison
genomes. On the shared server set:

```text
FASTQ_SCREEN_CONFIG="~micgdu/GenomicData/fastq_screen_db/fastq_screen.conf"
FASTQ_SCREEN_MISSING_ACTION="error"
```

The production `error` policy prevents a requested contamination/species check
from being silently skipped.

## Statistics and biological design

The default model and significance controls are:

```text
MIN_REPLICATES_PER_CONDITION=2
DESIGN="~ condition"
PAIRING_MODE="auto"
PAIRING_COLUMN="subject"
INCOMPLETE_PAIR_ACTION="error"
FDR=0.05
MIN_ABS_DELTA_PAU=0.10
```

`CONDITION_ORDER` controls deterministic pairwise contrast order; leave it
empty to use first appearance in the samplesheet. The workflow resolves and
records the actual paired or unpaired formula for every contrast before R is
run. Do not encode technical libraries or lanes as independent biological
replicates.

The Mcell2019 parameter block defines the registered APA-A method. Changing its
fixed discovery, internal-priming or gene-extension values creates a modified
method and must be documented in the project report.

## Parallel work and memory

Start with hard server-wide ceilings:

```text
MAX_TOTAL_THREADS=48
MAX_TOTAL_MEMORY_GB=384
```

Then tune per-stage job counts and per-job threads. The resource planner checks
both ordinary pools and the combined peak of concurrently running downstream
branches. For example, `STAR_PARALLEL_JOBS=4` with `STAR_THREADS=12` has a
48-thread ceiling before memory limits are considered.

```text
PARALLEL_DOWNSTREAM_MODULES=true
DOWNSTREAM_MODULE_PARALLEL_JOBS=3
STAR_PARALLEL_JOBS=4
STAR_THREADS=12
DGE_CONTRAST_PARALLEL_JOBS=3
APA_A_CONTRAST_PARALLEL_JOBS=4
APA_A2_CONTRAST_PARALLEL_JOBS=4
APA_B_CONTRAST_PARALLEL_JOBS=4
ENRICHMENT_PARALLEL_JOBS=6
TRACK_PARALLEL_JOBS=8
```

Always inspect `00_metadata/resource_plan.tsv` after validation. Reduce parallel
job counts before reducing the thread count required by an individual tool.
The [workflow-steps page](02_workflow_steps.md) explains which downstream
branches can overlap.

## Tracks and UCSC descriptors

Track families and normalization types are independently configurable. The
default publishes all-read, exact-end, filtered/rejected-end and active-PAS
families, including CPM and the applicable DESeq2/robust-CPM final signals.

```text
GENERATE_BIGWIGS=true
RETAIN_BEDGRAPH=false
UCSC_BIGDATA_URL_PREFIX=""
UCSC_NEGATE_MINUS_TRACKS=true
```

Leave the URL prefix empty when the public HTTPS address is not yet known. The
report still writes validated one-line descriptors using relative BigWig names;
the public prefix can be supplied in a later project configuration. Never put
Markdown link syntax in `UCSC_BIGDATA_URL_PREFIX`. See [tracks and
outputs](tracks_and_outputs.md).

## Cleanup

Success-only cleanup is on by default:

```text
CLEANUP_INTERMEDIATES=true
KEEP_TRIMMED_FASTQ=false
KEEP_LANE_BAMS=false
KEEP_TRACK_STRAND_BAMS=false
KEEP_BEDGRAPHS=false
```

Cleanup runs only after required receipts and final deliverables validate. It
does not remove final BAMs, count tables, BigWigs, reports or provenance. Enable
retention flags before a diagnostic run if those intermediates are needed.

## Validate every edited configuration

```bash
rna-ends2tracks --config /absolute/path/to/config.conf --stop-after validate
rna-ends2tracks --config /absolute/path/to/config.conf --dry-run
```

Review `validated_samples.tsv`, `contrasts.tsv`, `resource_plan.tsv` and
`warnings.tsv` before processing reads. Configuration accepted in one project
must not be assumed valid for another sample set, protocol or assembly.
