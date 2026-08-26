# Implementation status

## Alpha.6 source implementation

Implemented on `feature/workflow-revision` from the alpha.5 baseline:

- restricted `config.conf` plus one lane-level samplesheet;
- GRCh38/GRCm39 mixed-project separation and within-genome contrasts;
- bounded trimming, STAR, merge, extraction, track and contrast pools;
- C0-C5 count universes and invariant audits;
- exact Mcell2019 two-round discovery with documented tie/coordinate rules;
- installed-atlas rescue, atlas builder and provenance validation;
- C4 primary DESeq2 DGE, C5 diagnostics and per-contrast paired/unpaired models;
- DEXSeq APA shifts and intragenic candidate PCPA outputs;
- five track families with raw/CPM and final-signal DESeq2/robust CPM;
- locks, receipts, deterministic timing/index merges, reports and success-only cleanup;
- side-by-side shared-server installer and user/admin documentation.

## Validation still required before server promotion

Local source validation passes Python 3.11 compilation, 40 unit/portable-dry-run tests, Ruff static checks, parsing of all four R modules with R 4.4, Bash parsing, wheel construction and packaged-asset inspection. The laptop does not have the Linux bioinformatics executables or Bioconductor packages, and Windows has no compatible `pysam` wheel for the portable test interpreter; the tagged candidate must therefore also pass GitHub CI and the complete installed Linux environment tests. Before promoting alpha.6, build/install both versioned PAS atlases, run metadata/dry-run canaries for GRCh38 and GRCm39, repeat small real-read preprocessing/exact-end canaries, run a synthetic C0-C4 truth set and verify representative PAS/PCPA loci.

APA-B remains disabled until its independent pinned engine/model pilot is accepted. Alpha.6 must remain a prerelease until the above integration gates pass.
