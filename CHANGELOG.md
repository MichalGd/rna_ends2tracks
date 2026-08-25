# Changelog

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
