# C0-C5 data stages and count universes

The `C` labels are workflow-specific shorthand for **count universes**. They are not standard QuantSeq or RNA-seq terminology. The number identifies a processing level; `S` means a separated uncertain side branch and `R` means a rejected side branch. C5 is an independent diagnostic count matrix rather than the next transformation of C4.

| Label | Plain-language name | What it contains | Main use |
|---|---|---|---|
| **C0** | eligible aligned reads | Mapped, primary, uniquely assigned (`NH=1`) alignments. Duplicate flags are retained because the supported assay has no UMIs. | Final sample BAMs, all-read coverage tracks, and input to end extraction |
| **C1** | exact transcript-end counts | One strand-aware single-nucleotide transcript 3′ coordinate for each eligible C0 alignment whose end-defining CIGAR side is not clipped | Exact-end tracks and internal-priming evaluation |
| **C1S** | uncertain clipped-end counts | Eligible C0 alignments with soft/hard clipping at the side needed to define the exact 3′ coordinate | Separate QC; excluded from exact PAS calling |
| **C2** | filtered exact-end counts | C1 ends that pass the internal-priming mask, plus masked ends rescued by an assembly-matched annotated transcript end or PAS atlas | Condition-blind active-PAS discovery and filtered-end tracks |
| **C2R** | internal-priming rejects | C1 ends rejected by the A/T-rich internal-priming rule and not rescued | Diagnostic reject tracks; excluded from DGE and APA |
| **C3** | active-PAS counts | C2 ends counted at most once in the final non-overlapping, project-wide active-PAS intervals | PAS usage matrix, APA tests, and active-PAS tracks |
| **C4** | active-PAS gene counts | C3 counts summed per gene over uniquely assigned PAS | Primary raw-integer gene-expression matrix and normalization factors |
| **C5** | conventional exon-count diagnostic | Reverse-stranded featureCounts exon-overlap counts from C0 BAMs | Diagnostic comparison with C4; never substituted for C4 |

## Required count relationships

For every sample, the workflow audits:

```text
C0 = C1 + C1S
C1 = C2 + C2R
sum(C3) <= sum(C2)
C4 = gene-wise sums of uniquely assigned C3
```

APA-B may reuse receipt-validated C1 and C1S as an execution optimization. It adds the two tables back together, verifies the recorded `C0 = C1 + C1S` identity and source receipt, and therefore reconstructs the same raw endpoint universe that its BAM reader would produce. The exact-end audit also records the duplicate-flagged C0 count; reads remain retained because this is a no-UMI workflow. This does not reuse C2, active PAS, gene assignments, or APA-A statistics; PolyAseqTrap/DeepIP discovery remains an independent method.

C3 can be smaller than C2 because filtered ends outside the discovered active-PAS intervals are not assigned to C3. Ambiguous, antisense, and unassigned C3 sites can remain in catalogs and tracks but are excluded from C4.

C5 has no arithmetic equality with C4. The two matrices answer different questions: C4 counts reads at accepted active PAS, whereas C5 counts conventional annotated-exon overlaps. A large C4/C5 difference is reported for review, not automatically corrected by replacing one matrix with the other.

## Short example

Suppose a sample has 100 eligible C0 alignments. Five have an end-defining clip, so C1 contains 95 exact ends and C1S contains 5 uncertain ends. If 10 exact ends fail the internal-priming rule, C2 contains 85 and C2R contains 10. If 80 C2 ends fall inside active-PAS intervals, C3 totals 80. Only the uniquely gene-assigned portion of those 80 is summed into C4.

See [methods](methods.md) for coordinate and filtering rules and [tracks and outputs](tracks_and_outputs.md) for the track generated from each universe.
