# FastQ Screen contamination and species QC

FastQ Screen is an optional, non-filtering quality-control step. It screens a configurable subset of every raw lane and every available mate against the genomes and contaminant references defined by the server's FastQ Screen configuration. Reads are not removed and the screening result never changes alignment, DGE, or APA counts.

## Server requirement

The workflow environment provides `fastq_screen` and Bowtie2, but the large reference indexes are deliberately site-managed and are not bundled with a release. An administrator should create one readable FastQ Screen configuration containing the expected species plus appropriate controls/contaminants, for example human, mouse, PhiX, vectors/adapters, rRNA, bacterial references, and any organisms routinely processed by the facility.

Set the path in `config.conf`:

```text
RUN_FASTQ_SCREEN=true
FASTQ_SCREEN_CONFIG="/shared/references/fastq_screen/fastq_screen.conf"
FASTQ_SCREEN_MISSING_ACTION="warn"
FASTQ_SCREEN_SUBSET=200000
FASTQ_SCREEN_THREADS=4
FASTQ_SCREEN_PARALLEL_JOBS=4
FASTQ_SCREEN_MEMORY_GB=4
```

`FASTQ_SCREEN_MISSING_ACTION="warn"` preserves portability: the workflow records `SKIPPED_MISSING_CONFIG` and continues. Use `"error"` on a production server where contamination screening is mandatory. A skipped check is never reported as a pass.

For paired-end libraries, R1 and R2 are both screened and reported in separate mate subdirectories. This also prevents output-name collisions when mate files happen to have the same basename in different source directories. FastQ Screen treats the files as independent QC inputs; it does not alter pairing or produce the alignment input.

## Outputs and interpretation

Outputs are retained under `01_qc/fastq_screen/<sample>.<technical-replicate>.<lane>/`. `fastq_screen_summary.tsv` records the layout, screened mates, status, output directory, and exact configuration path. `fastq_screen_metrics.tsv` consolidates every mate/database row, including processed-read and mapping-percentage columns. Native text/PNG/HTML reports are collected by MultiQC and linked in the final report. `PASS` means technical execution succeeded; it is not an automatic declaration that the sample is uncontaminated.

Interpretation depends on the databases in the site configuration. A mouse library should be dominated by mouse-compatible sequence and a human library by human-compatible sequence; cross-mapping between closely related genomes is expected and is not itself proof of contamination. Review unexpected organism, vector, PhiX, or microbial signal together with FastQC, STAR mapping, and sample provenance.
