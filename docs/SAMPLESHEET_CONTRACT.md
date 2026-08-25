# Samplesheet metadata contract

Each row describes one single-end R1 FASTQ from one sequencing lane. A project
contains one genome assembly and may contain many conditions, biological samples,
technical library preparations, and lanes.

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
`mm39`. Aliases are normalized to the canonical assembly name. All rows in one
project must resolve to the same assembly, which must exactly match `assembly` in
the selected reference manifest. A mismatch stops validation before FASTQ
processing.

The manifest—not the samplesheet alone—selects the actual FASTA, GTF, STAR index,
chromosome sizes, and optional PAS atlas. This prevents a short genome label from
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

## Migration from alpha.4

Alpha.5 makes `description`, `genome`, and `technical_replicate_id` required
columns. For an existing sample with one library and one lane, retain its current
`sample_id`, `biological_replicate_id`, and `lane_id`, set
`technical_replicate_id` to a stable value such as `T01`, and set `genome` to the
assembly in the reference manifest.

Do not change workflow release within an active project. Alpha.5 adds technical
library IDs to intermediate filenames and read groups, so an alpha.4 project must
finish with its versioned alpha.4 launcher. Start a new output directory when
adopting the alpha.5 contract.
