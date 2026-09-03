# APA-A2 corrected analysis

APA-A2 is a separate, independently switchable correction of the legacy APA-A
statistical analysis. It does not replace, rewrite, or suppress APA-A. New
projects run APA-A, APA-A2, and validated APA-B by default.

## What is independent

APA-A2 has its own workflow stage (`apa_a2`), DEXSeq invocation, logs,
per-contrast receipts, module receipt, result index, enrichment jobs, output
directory, and report columns. Disable it with `RUN_APA_A2=false` without
changing APA-A or APA-B.

APA-A2 and APA-A deliberately test the same condition-blind C3 active-PAS
catalog. This makes APA-A2 a corrected reanalysis of the same biological site
universe, not a third site-discovery algorithm. APA-B remains the method with
an independently discovered PolyAseqTrap/DeepIP catalog.

## Corrected effect calculation

DEXSeq is fitted independently from raw C3 counts using the resolved paired or
unpaired design. APA-A2 then calculates effects directly from the raw count
matrix:

1. within every gene and sample, PAS usage (PAU) is the PAS count divided by
   that gene's total selected-PAS count;
2. a gene/sample with zero total counts is unavailable (`NA`), never assigned
   artificial zero PAU;
3. an unpaired effect is mean numerator PAU minus mean denominator PAU;
4. a paired effect is calculated within every complete subject pair and then
   averaged with equal weight across pairs;
5. all effect vectors are explicitly aligned by `pas_id` before publication.

The primary APA-A2 site definition requires both `padj <= FDR` and
`abs(delta_PAU) >= MIN_ABS_DELTA_PAU`. A primary gene requires a DEXSeq
gene-level adjusted p-value at or below FDR and at least one primary site.
Proximal/distal direction is derived from the PAU-weighted transcript-coordinate
shift over the complete selected site set; the sign is corrected for transcript
strand.

## Outputs

Each genome is written under `06b_apa_a2_corrected/<genome>/`:

- `dexseq_a2/<contrast>.apa_a2_sites.tsv`: independent DEXSeq statistics,
  condition PAUs, corrected delta-PAU, and significant/primary flags;
- `dexseq_a2/<contrast>.apa_a2_genes.tsv`: gene q-values, effect-qualified
  site counts, weighted shift and proximal/distal class;
- `dexseq_a2/<contrast>.apa_a2_pair_deltas.tsv`: auditable within-subject
  effects for paired designs (header-only for unpaired designs);
- `dexseq_a2/<contrast>.apa_a2_audit.tsv`: PAU-sum, effect-conservation,
  missing-observation, and design audit;
- `candidate_pcpa.tsv`: effect-qualified candidate premature
  cleavage/polyadenylation events;
- `enrichment/`: APA-A2 GO, Reactome, Hallmark, and KEGG ORA/GSEA outputs.

`10_reports/alternative_polyadenylation_summary.tsv`,
`top_apa_gene_events.tsv`, and the Markdown/HTML report include APA-A2. The
existing `08_apa_comparison/` contract remains the independent-catalog APA-A
versus APA-B comparison and is not silently redefined.

## Interpretation

Legacy APA-A remains available for continuity and method comparison. APA-A2
primary calls should be preferred when an effect-size-qualified interpretation
of the Mcell2019 C3 universe is required. Large discrepancies between APA-A and
APA-A2 are informative: they may reflect the legacy normalized-count PAU
calculation or its two-site direction heuristic, not a change in the underlying
C3 catalog or DEXSeq hypothesis test.
