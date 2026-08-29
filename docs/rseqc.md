# RSeQC quality control

RSeQC is a default, receipt-backed post-alignment stage. It analyzes each final C0 BAM after STAR/sample merging and before the end-derived analyses. It does not alter reads, PAS calls, gene counts, or statistical models.

The implementation follows the [RSeQC documentation](https://rseqc.sourceforge.net/) and emits the native filenames recognized by the [MultiQC RSeQC module](https://docs.seqera.io/multiqc/modules/rseqc/).

## Analyses

| Tool | Purpose | Main interpretation for QuantSeq REV |
|---|---|---|
| `infer_experiment.py` | Independent annotation-aware strandedness estimate | The reverse-compatible fraction should dominate for the validated REV profile |
| `read_distribution.py` | Distribution among annotated exons, introns and flanking regions | Diagnostic context for library composition; not a pass/fail measure by itself |
| `geneBody_coverage.py` | Relative coverage across transcript percentiles from 5-prime to 3-prime | Strong 3-prime enrichment is expected for a 3-prime counting assay |

The workflow records the first and last 20 gene-body percentiles as 5-prime and 3-prime fractions and reports their ratio. `quantseq_three_prime_enriched=true` means the last 20 percentiles carry more signal than the first 20. It is a transparent descriptive flag, not a universal biological threshold.

## Annotation contract

RSeQC requires BED12. Set `HG38_RSEQC_BED` and/or `MM39_RSEQC_BED` to audited assembly-matched BED12 files, or leave them empty. Empty values cause the workflow to create a deterministic transcript BED12 file from the exact GTF already selected for that genome. This prevents accidental mixing of assemblies or annotation releases. The automatically generated BED contains all GTF transcripts; an explicitly configured, validated subset such as a housekeeping-gene BED12 remains possible when that is scientifically preferred.

## Configuration

```bash
RUN_RSEQC=true
RSEQC_INFER_EXPERIMENT=true
RSEQC_READ_DISTRIBUTION=true
RSEQC_GENE_BODY_COVERAGE=true
RSEQC_MULTIQC=true
RSEQC_SAMPLE_READS=200000
RSEQC_MIN_TRANSCRIPT_LENGTH=100
RSEQC_PARALLEL_JOBS=6
RSEQC_MEMORY_GB=4
HG38_RSEQC_BED=""
MM39_RSEQC_BED=""
```

At least one analysis must remain enabled when `RUN_RSEQC=true`. `RSEQC_PARALLEL_JOBS` controls independent per-sample jobs and is checked against the global thread/memory ceilings.

## Outputs

```text
01_qc/rseqc/
├── references/                       generated BED12 files, when needed
├── infer_experiment/                 per-sample strandedness reports
├── read_distribution/                per-sample annotation-distribution reports
├── gene_body/                        native per-sample RSeQC tables/plots
├── gene_body_coverage.tsv            normalized 100-percentile matrix
├── gene_body_coverage.svg            combined cross-sample plot
├── rseqc_summary.tsv                 concise per-sample metrics
├── multiqc/multiqc_report.html       dedicated RSeQC MultiQC report
└── run_receipt.json                  resume/integrity receipt
```

The final `10_reports/report.html` embeds the combined plot, presents the per-sample summary, and links the detailed RSeQC/MultiQC artifacts.

## Practical interpretation

A successful QuantSeq REV experiment should combine: strong STAR mapping, the expected reverse-compatible orientation, a clear 3-prime gene-body bias, biological replication, and plausible count/PAS summaries. RSeQC cannot by itself prove that a library is biologically successful, and its gene-body curve must not be judged against the near-uniform expectation used for conventional full-length RNA-seq.
