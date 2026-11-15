# Gene expression in alpha13

Alpha13 defaults to `extended_exons_v1`, independent of PAS discovery. Existing PAS, APA, alignment and track algorithms retain their definitions. Legacy C4 differential-expression calculations are preserved and opt-in; alpha13.post1 changes only its distance-plot presentation.

In alpha13.post1, both extended-exon and legacy C4 sample-distance plots use complete-linkage hierarchical clustering of the original Euclidean VST distances. The same clustered sample order appears left to right and top to bottom, with a top-left diagonal. Labels no longer use alphabetical order, and figure dimensions grow with sample count. Numeric distance tables and differential-expression calculations are unchanged. The bundled `scripts/replot_vst_distances.R` accepts a saved distance TSV and output prefix to regenerate just the PNG/PDF without alignment, counting or model fitting.

The default method settings are:

```ini
DGE_METHODS=extended_exons_v1
DGE_PRIMARY_METHOD=extended_exons_v1
# To include legacy DGE: DGE_METHODS=extended_exons_v1,legacy_c4
# For legacy-only, set both keys to legacy_c4.
RUN_DGE_STRINGENT=true
RUN_DGE_RELAXED=true
DGE_STRINGENT_LFC_THRESHOLD=1.0
FDR=0.05
GEX_EXTENSION_MAX_NT=6000
GEX_TSS_BUFFER_NT=2000
GEX_MIN_COUNT=5
GEX_MIN_TOTAL_COUNT=15
```

## Annotation and counts

The source assembly-matched GENCODE GTF stays unchanged. A dedicated SAF contains the union of original exons per gene plus a contiguous extension beyond the outermost annotated exon in transcriptional direction. Shorter isoforms are represented by their existing exons; their ends are not extended into the same gene's introns. Full gene IDs, including versions, are counting keys.

The default cap is 6,000 nt. Another same-strand gene's body or upstream TSS buffer stops the extension, including nested/overlapping genes and all alternative annotated TSSs. If no transcript TSS is present, the gene boundary is used and audited. A buffer already covering the proposed extension adds zero bases; original exons are never shortened. Chromosome boundaries and remaining same-strand extension conflicts also stop extension. Opposite-strand overlaps are retained under stranded counting and recorded in the audit. Internal coordinates are zero-based half-open; exported SAF is one-based inclusive.

Extended biotypes are protein_coding, lncRNA, lincRNA, antisense, processed_transcript, sense_intronic, sense_overlapping, 3prime_overlapping_ncRNA, macro_lncRNA, and IG/TR C/D/J/V gene categories listed in `gex_annotation.EXTEND_BIOTYPES`. Original exons of other biotypes are retained without extension; every biotype acts as a neighbor barrier. The manifest records the complete allowlist and hashes.

Counting uses C0 mapped-primary NH=1 sample BAMs and `featureCounts -F SAF -s 2 --primary`. Paired-end libraries additionally use `-p --countReadPairs`; biological pairing is a separate design choice. Duplicate flags are retained. No fractional/multimapper/multi-gene counting is enabled. Overlapping features of one gene count once, and ambiguous overlaps across compatible-strand genes are excluded. The complete matrix remains integer and unfiltered. [featureCounts documentation](https://subread.sourceforge.net/featureCounts.html)

The 6-kb/2-kb distances are a project policy, not prescribed Lexogen defaults. Exon-union plus extension is broader than terminal-exon-only counting. It can recover incomplete 3-prime annotation, but also internal priming or readthrough; low-count filtering does not diagnose these artifacts. New DGE bypasses PAS/end-clipping/internal-priming eligibility gates. PAS analyses keep those gates. Setting extension to zero provides a conventional exon-union comparator. [Lexogen QuantSeq REV V2 manual, page 30](https://www.lexogen.com/wp-content/uploads/2026/01/225UG675V0102_QuantSeq-REV-V2-3-mRNA-Seq_2026-01-19.pdf#page=30)

## Filtering and two statistical tests

Each contrast filters its own raw count subset with `edgeR::filterByExpr(min.count=5, min.total.count=15)`, using prefilter assigned-count library sizes and the resolved design matrix. In simple small-group designs, the CPM rule requires expression equivalent to five counts at the median library depth in N samples, where N is the smaller group size. Support is counted across both groups, retaining genes expressed predominantly in one condition. Complex designs use design leverage; groups above ten use edgeR's `large.n=10, min.prop=0.7` adjustment. Use `GEX_MIN_COUNT=10` for an explicit sensitivity rerun. [edgeR manual](https://bioconductor.org/packages/release/bioc/manuals/edgeR/man/edgeR.pdf)

One DESeq2 model is fitted per contrast with poscounts size factors. Raw integer counts enter the model. Automatic independent filtering is disabled after the explicit prefilter, and both result tables share the finite valid-test universe; Cook's/outlier safeguards remain enabled. BH adjustment is performed separately for each hypothesis over that common universe:

| Mode | Hypothesis rejected at adjusted p <= FDR | Default |
|---|---|---|
| `stringent` | absolute log2FC <= 1, using `lfcThreshold=1, altHypothesis="greaterAbs"` | On |
| `relaxed_fdr` | log2FC = 0 | On |

Stringent means evidence that the effect exceeds twofold, not ordinary significance followed by a fold-change cutoff. Both modes use unshrunken MLE fold changes in tables/plots. No global FDR across modes or contrasts is claimed. Large DE fractions can reflect cell-type differences or normalization/design problems; adding an effect threshold does not validate the model. [DESeq2 vignette](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html)

## ORA, GSEA and report

The report orders stringent DE/ORA, relaxed DE/ORA, then shared GSEA. ORA uses each mode's significant genes with separate higher-numerator, higher-denominator and nondirectional `any_dge` queries. Relaxed ORA has no additional fold-change cutoff. Legacy C4 retains its original foreground rule.

Shared GSEA runs once from all common valid genes ranked by the ordinary zero-null Wald statistic, even if one or both DE modes are disabled. It never uses a significant-only subset or the minimum-effect statistic. Positive NES favors the numerator and negative NES favors the denominator. Summaries select terms separately by collection and direction; full tables remain available. [fgsea tutorial](https://bioconductor.org/packages/release/bioc/vignettes/fgsea/inst/doc/fgsea-tutorial.html)

ORA BH correction is per query/database/contrast/mode; GSEA BH correction is per database/contrast. ORA's background includes eligible valid genes absent from the selected databases, consistent with the previous denominator convention; mapping metrics and per-gene coverage make this explicit. Database coverage and redundant terms require interpretation. Enrichment alone is not evidence of biochemical pathway activation. MSigDB source/version, mouse orthology mapping and fixed seed are recorded.

```
05_gene_expression/<genome>/extended_exons_v1/
  annotation/                 SAF, BED, gene metadata, extension audit, hashes
  raw_counts.tsv              complete gene matrix
  featureCounts.tsv.summary   assignment categories
  qc/                         global size factors, VST PCA and sample distances
  result_index.tsv            explicit method/mode identity
  contrasts/<contrast>/
    model.rds                 shared fitted model with input/version signature
    eligibility.tsv           every input gene and eligibility reason
    filter_parameters.tsv     depths and per-sample thresholds
    test_universe.tsv          valid/invalid tests after model fitting
    ranking.tsv               ordinary Wald results for shared GSEA
    normalized_counts.tsv     eligible genes, per-contrast normalization
    stringent/                DE tables, MA/volcano, enrichment/ (ORA)
    relaxed_fdr/              DE tables, MA/volcano, enrichment/ (ORA)
    shared_gsea/enrichment/    ranked GSEA
```

Legacy paths remain `C4_primary_deseq2/` and `enrichment/<contrast>/`. Disabled new modes are omitted from indexes/report; existing mode files are retained as historical files, not deleted. New-method report figures and linked tables are copied below `10_reports/assets/`, with checksums in the asset manifests. Download the whole `10_reports` folder. Existing generic links to large upstream/BAM/browser files still require those originals.

Counts, QC and per-contrast results have independent receipts. Changing the expression filter reuses counts and refits statistics. Changing only the two DE switches, FDR or effect threshold reuses the compatible fitted model. Annotation changes rebuild counts and downstream results. Old C4 outputs are never accepted as new-method counts. Alpha13 does not grant blanket cross-version receipt compatibility.

## Reuse final sample BAMs

Copy the old samplesheet, retaining sample IDs, library metadata, contrast direction/order and biological design. Historical FASTQ paths may remain in it; BAM mode does not require those files. Use a new output directory and the [BAM rerun example](../config/expression_bam.conf.example).

Alternatively, `rna-ends2tracks-prepare-expression --source-config /old/config.conf --new-project /dysk3/micgdu/Analysis/run11_quanSeq/mouse_gex_v1 --project-id mouse_gex_v1 --jobs 8` prepares the config, copied samplesheet and source metadata without launching work. It refuses a nonempty destination, preserves pairing/design/condition order, and checks available old validated contrast definitions for agreement. It applies the new expression defaults; review the prepared files before the later server run.

`INPUT_MODE=bam` reads `<ALIGNMENT_SOURCE_DIR>/<sample_id>/<sample_id>.bam` in place, with BAI/CSI. This first interface supports the new expression branch plus enrichment/report only: disable APA, RSeQC, tracks and cleanup. It validates file completeness, coordinate-sorted header, assembly dictionary, available sample read-group identity, index, a sample of 10,000 C0-policy records, and available merge-receipt size/hash/timestamp evidence. The policy check is sampled, not a full independent audit of every alignment; use original workflow C0 BAMs. Any missing BAM is an error, never a trigger to rerun FASTQs. External inputs are not copied, modified or cleanup outputs.

After installing alpha13 in a separate environment:

```bash
rna-ends2tracks --config /path/to/new/config.conf --stop-after validate
rna-ends2tracks --config /path/to/new/config.conf
```

All new count/statistical/report outputs follow `OUTPUT_DIR`; new expression/enrichment temporary files follow `TMP_DIR` (or the new output tree). Set both under `/dysk3`. No STAR/QC/trimming/PAS work is scheduled in BAM mode. The new expression tests can use eight contrast jobs; STAR thread settings have no effect during a BAM-only rerun.

Before production, use a short comparison with zero-extension counts, inspect recovered genes and ambiguity, confirm replicate clustering/normalization and selected marker loci, and verify actual Linux BAM input validation. There is no acceptance target for total tested genes or significant-gene percentage. Existing mouse/human outputs and their BAMs should remain unchanged.
