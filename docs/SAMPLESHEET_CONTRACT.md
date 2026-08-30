# Samplesheet metadata contract

Each row describes one sequencing lane. For `SE`, `fastq_r1` is required and
`fastq_r2` must be empty. For `PE`, both mate files are required and must be
distinct. A project uses one library layout consistently. It may contain
GRCh38, GRCm39, or both, and may contain many conditions, biological samples,
technical library preparations, and lanes.

Supported protocol values are `quantseq_rev_v1_se`, `quantseq_rev_v2_se`,
`quantseq_rev_v1_pe`, and `quantseq_rev_v2_pe`; the suffix must match
`library_layout`. `umi_present` must be false. In paired mode both mates are
used for alignment and conventional coverage, but QuantSeq REV R1 alone is the
end-defining mate for APA counting.
The default preprocessing removes the first 12 random-primer-derived bases
from PE R2 only (`PE_R2_TRIM_5P=12`); R1 is not force-trimmed.

The repository includes `config/samplesheet.example.csv` for SE and
`config/samplesheet.paired_end.example.csv` for PE. For a paired project also
set `LIBRARY_LAYOUT="paired_end"` and `LIBRARY_PROTOCOL="quantseq_rev_v2_pe"`
in `config.conf`.

## Replicate hierarchy

| Column | Meaning | Statistical role |
|---|---|---|
| `sample_id` | One biological analysis unit and final merged BAM | One sample column in DGE and APA count matrices |
| `biological_replicate_id` | Independent biological specimen, culture, animal, or preparation | Counts toward the minimum biological replication for a condition |
| `technical_replicate_id` | Independently prepared library from the same biological unit | Merged into `sample_id`; never counted as an independent replicate |
| `lane_id` | Sequencing lane for one technical library | Merged into `sample_id`; never counted as an independent replicate |
| `subject` | Biological unit shared across conditions in a genuinely matched design | Activates a paired formula only when matching is complete for a contrast |

All rows sharing `sample_id` must have identical biological metadata, including
description, genome, biological replicate, condition, batch, subject, protocol,
read length, kit, and UMI status. They may differ in `technical_replicate_id`,
`lane_id`, and FASTQ path.

The tuple `(sample_id, technical_replicate_id, lane_id)` must be unique. This
allows one technical library to span several lanes and also allows separate
technical libraries to use the same lane label.

## Genome selection

The samplesheet accepts `GRCh38` or its alias `hg38`, and `GRCm39` or its alias
`mm39`. Aliases are normalized to the canonical assembly name. Each row selects its
assembly-specific reference group from `config.conf`. Processing, condition-blind PAS
discovery and contrasts remain genome-specific; cross-genome contrasts are forbidden.
A missing or mismatched asset stops validation before FASTQ processing.

`config.conf`—not the samplesheet alone—selects the actual FASTA, GTF, STAR index,
chromosome sizes, and PAS atlas. This prevents a short genome label from
silently selecting an unversioned or incompatible reference.

## Examples

One biological sample with two technical libraries, where the first library was
sequenced on two lanes, uses three rows:

```csv
sample_id,description,genome,biological_replicate_id,technical_replicate_id,lane_id
CTRL_01,Control replicate 1,GRCh38,CTRL_01,T01,L001
CTRL_01,Control replicate 1,GRCh38,CTRL_01,T01,L002
CTRL_01,Control replicate 1,GRCh38,CTRL_01,T02,L001
```

Those rows generate three lane BAMs with distinct read groups and one merged
`CTRL_01.bam`. Downstream DGE and APA receive one `CTRL_01` count column.

For paired control/treatment samples from the same animal or culture lineage, use
different `sample_id` and `biological_replicate_id` values but repeat `subject`:

```csv
sample_id,biological_replicate_id,condition,subject
CTRL_01,CTRL_01,control,SUBJECT_01
TRT_01,TRT_01,treatment,SUBJECT_01
```

Technical replication must never be encoded by duplicating a biological sample
under new `sample_id` values. Doing so would pseudoreplicate the statistical
analysis.

## Migration to alpha.6

Alpha.6 retains the alpha.5 columns, uses `config.conf` for normal project settings,
and permits mixed-genome metadata. For an existing sample with one library and one lane, retain its current
`sample_id`, `biological_replicate_id`, and `lane_id`, set
`technical_replicate_id` to a stable value such as `T01`, and set `genome` to the
assembly configured in `config.conf`.

Do not change workflow release within an active project. Finish older projects with
their versioned launcher and start a new output directory when adopting alpha.6.
