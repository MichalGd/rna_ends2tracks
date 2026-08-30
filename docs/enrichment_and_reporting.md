# Statistical plots, enrichment, provenance, and APA-B interpretation

## What the workflow produces

The primary C4 DESeq2 analysis writes variance-stabilized PCA coordinates, a PCA figure, a sample-distance matrix and heatmap. Every contrast writes its complete DESeq2 table, MA plot, and volcano plot in PDF and PNG. PNG files are embedded in `10_reports/report.html`; PDF and source TSV files remain linked for publication-quality review.

Gene-set interpretation is a separate workflow stage after DGE and APA statistics. It never changes a DGE or APA result. The stage prepares one row per tested gene, then runs:

- over-representation analysis (ORA) on configured significant foregrounds;
- exploratory ranked GSEA using the complete testable-gene background;
- GO Biological Process, Molecular Function, and Cellular Component collections;
- Reactome and Hallmark collections;
- KEGG MEDICUS and KEGG LEGACY collections from the installed, version-recorded MSigDB snapshot.

KEGG enrichment uses the gene sets exposed by the pinned `msigdbr`/MSigDB installation; it does not query the live KEGG API or redistribute KEGG pathway files. KEGG MEDICUS is preferred, while KEGG LEGACY is retained as a separately labelled database. The exact MSigDB and package versions are written to every job's `provenance.tsv`.

## Analysis-specific universes

The analyses are intentionally independent:

| Analysis | Tested background | ORA foregrounds | GSEA ranking |
|---|---|---|---|
| DGE | genes tested by the contrast-specific C4 DESeq2 model | any DGE, upregulated, downregulated | DESeq2 Wald statistic; signed p-value fallback only when the statistic is unavailable |
| APA-A | genes tested by DEXSeq with at least two active PAS | any APA, distal shift, proximal shift | signed gene-level DEXSeq shift score |
| APA-A PCPA | same APA-A testable-gene universe | any candidate PCPA, increased PCPA, decreased PCPA | signed APA shift score |
| APA-B | genes tested by the independent DRIMSeq/stageR model | the corresponding APA and PCPA directions | signed gene-level APA-B shift score |

Mouse collections are obtained through `msigdbr`'s recorded human-to-mouse orthology mapping. Each job writes `mapping_audit.tsv` and `provenance.tsv`, including package and MSigDB version information. Mapping loss, background size, and foreground size are visible in the final report. An empty enrichment result is valid and produces header-only tables plus an explicit no-significant-term plot.

## Configuration

```text
RUN_DGE_ENRICHMENT=true
RUN_APA_ENRICHMENT=true
ENRICHMENT_ORA=true
ENRICHMENT_GSEA=true
ENRICHMENT_GO=true
ENRICHMENT_REACTOME=true
ENRICHMENT_HALLMARKS=true
ENRICHMENT_KEGG=true
ENRICHMENT_RICH_PLOTS=true
ENRICHMENT_NETWORK_MAX_TERMS=8
ENRICHMENT_NETWORK_MAX_GENES=50
ENRICHMENT_PADJ=0.05
ENRICHMENT_DGE_MIN_ABS_LFC=1.0
ENRICHMENT_APA_MIN_ABS_DELTA_PAU=0.10
ENRICHMENT_MIN_GENESET_SIZE=10
ENRICHMENT_MAX_GENESET_SIZE=500
ENRICHMENT_PARALLEL_JOBS=6
```

`FDR` controls selection of DGE/APA effects. `ENRICHMENT_PADJ` controls pathway significance. These thresholds answer different questions and are recorded separately.

## Outputs

Every enrichment job contains `prepared_gene_table.tsv`, `ora.tsv`, `gsea.tsv`, `mapping_audit.tsv`, `enrichment.pdf`, `enrichment.png`, `provenance.tsv`, `plot_index.tsv`, and a receipt. When `ENRICHMENT_RICH_PLOTS=true`, the `plots/` subdirectory contains database-specific ORA and GSEA dotplots, barplots, and gene-to-concept network plots in both PNG and PDF. Networks are deliberately capped by `ENRICHMENT_NETWORK_MAX_TERMS` and `ENRICHMENT_NETWORK_MAX_GENES` to remain legible. The global machine-readable index is `10_reports/enrichment_summary/enrichment_index.tsv`.

The scientific report publishes both overview and drill-down tables:

- `differential_gene_expression_summary.tsv`: tested, significant, up- and downregulated genes per contrast;
- `top_differential_genes.tsv`: up to 25 FDR-significant genes per contrast with effect sizes and adjusted p-values;
- `alternative_polyadenylation_summary.tsv`: APA-A, APA-B, PCPA and method-concordance counts per contrast;
- `top_apa_gene_events.tsv`: significant gene-level events from each independent APA method, including shift direction and PCPA status;
- `top_enrichment_terms.tsv`: the strongest significant ORA and GSEA terms for every method/contrast job.

These are summaries, not replacements for the complete contrast tables referenced by each module's `result_index.tsv`. The HTML report embeds the same summaries and explains thresholds and interpretation limits.

The final report also creates `10_reports/provenance_dashboard/`:

- `provenance_summary.tsv`: workflow, input, reference, PAS-atlas, and APA-B identity;
- `receipt_inventory.tsv`: every upstream receipt and its signature/status;
- `environment_packages.tsv`: complete Conda package inventory when available;
- `software_versions.tsv`: executable paths and version probes;
- `output_manifest.tsv`: every output, size, modification time, and SHA-256 for small files;
- `dashboard.json`: compact machine-readable run overview.

Large files use size and modification time in the manifest to avoid re-reading hundreds of gigabytes during report generation. Their stage receipts apply the same audited policy.

## APA-B interpretation gate

APA-B is not inferred from APA-A. It can be interpreted only when all of the following are true:

1. `RUN_APA_B=true` and `APA_B_PILOT_ACCEPTED=true`;
2. `APA_B_VALIDATION_MANIFEST` is an accepted schema-v1 manifest;
3. engine commit, model checksum, and environment checksum are pinned;
4. synthetic strand/coordinate tests and real QuantSeq REV canaries pass for every requested assembly;
5. the adapter's run-specific `engine_provenance.json` matches the accepted manifest;
6. UMI processing and coordinate deduplication are both explicitly false.

Otherwise the report says `DISABLED_NOT_VALIDATED` and does not fabricate results. SE and PE are separate validation claims: a paired project requires its exact `quantseq_rev_*_pe` profile in `library_protocols` and a real paired-layout canary in the accepted manifest. When enabled, APA-A and APA-B keep separate catalogs, gene summaries, enrichment results, and PCPA calls. Only the dedicated comparison output reports positional and direction concordance.

This is a scope check, not a statement that APA-B is unfinished. The audited `biolserv` GRCh38/GRCm39 QuantSeq REV V2 single-end deployment satisfies the gate when its combined accepted manifest is configured. Paired-end APA-B does not inherit that acceptance automatically.
