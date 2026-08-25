# Changelog

## 0.1.0-alpha.4 - 2026-08-25

- Added configurable contrast-specific pairing with `none`, `auto`, and `required` modes.
- Resolve complete matched-subject contrasts to `~ subject + condition` while retaining the default unpaired formula for contrasts with disjoint subject sets.
- Reject partially matched subjects by default and validate every resolved pair-specific design for completeness and full rank before compute.
- Reject one-level design terms before invoking R, avoiding late `model.matrix` factor errors.
- Fit DESeq2 independently per contrast and apply the same resolved formula to DEXSeq and DRIMSeq/stageR.
- Record pairing status, number of pairs, design mode, and resolved formula in metadata, statistical indexes, reports, and receipts.
- Isolated preprocessing signatures from statistical-only configuration changes so later statistical edits do not force unnecessary preprocessing reruns.

## 0.1.0-alpha.3 - 2026-08-25

- Removed the incompatible Anaconda `defaults` channel and explicitly disabled implicit default channels.
- Updated the samtools pin from 1.20 to 1.21 to satisfy the BBMap 39.13 runtime contract.
- Added a Linux CI job that creates the complete Conda environment with strict Bioconda/conda-forge channel priority and verifies required Python, R and command-line tools.
- Restricted BBDuk quality trimming to the right end of REV Read 1 so preprocessing cannot move its APA-defining 5-prime coordinate.

The corrected specification was first dry-run solved with Mamba 2.8.0 on the target `biolserv` server: 386 packages and approximately 450 MB of downloads.

## 0.1.0-alpha.2 - 2026-08-25

- Fixed the non-dry-run preprocessing order so the QuantSeq REV orientation check reads STAR gene-count output only after STAR finishes.
- Added regression coverage for a complete mocked single-end lane through alignment, orientation validation, merge and receipt handoff.
- Documented reuse of one audited `sjdbOverhang=149` index per species for read lengths up to 150 bp.

This remains an alpha canary. Linux tool integration, synthetic truth-set validation and the PolyAseqTrap REV pilot are still pending.

## 0.1.0-alpha.1 - 2026-08-25

Initial pre-release implementation:

- shared no-UMI Lexogen QuantSeq REV single-end preprocessing and alignment;
- human GRCh38 and mouse GRCm39 reference contracts;
- independent featureCounts/DESeq2 gene-expression analysis;
- direct-end APA-A with DEXSeq and candidate intragenic PCPA classification;
- pilot-gated PolyAseqTrap/DRIMSeq/stageR APA-B interface;
- cross-method comparison, strand tracks, receipts, schemas and server documentation.

This alpha has not completed Linux tool integration, public-data validation or the PolyAseqTrap REV pilot. It must not be described as production validated.
