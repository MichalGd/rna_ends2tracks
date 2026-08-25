# Implementation status

## Implemented

- Portable `rna-ends2tracks` CLI with independently callable modules and an `all` sequence.
- Strict lane-level samplesheet with explicit assembly, biological replicate, technical library, sequencing lane, description, no-UMI protocol, additive-design and reference validation.
- Deterministic all-pairs contrasts for conditions with at least two biological replicates; n=2 comparisons are labelled exploratory.
- Configurable contrast-specific pairing with complete-subject detection, mixed paired/unpaired projects, pair-specific rank validation, and resolved-formula provenance.
- Human GRCh38 and mouse GRCm39 reference-manifest profiles using the same analysis code.
- Lexogen QuantSeq REV SE raw-FASTQ QC, adapter/poly(A)/poly(T)/quality trimming, technical-library/lane read groups, STAR alignment, technical-library/lane merging, mapping QC and empirical orientation checks.
- Independent featureCounts and pair-specific DESeq2 gene-expression branch.
- APA-A direct-end extraction, audit, internal-priming evidence, condition-blind PAS clustering, annotation, raw counts, DEXSeq/Delta-PAU, and PCPA classification with a tested terminal comparator.
- Pilot-gated APA-B adapter contract, independent output validation, DRIMSeq/stageR, and branch-local PCPA classification.
- APA catalog/effect/PCPA comparison, transcript-strand BigWigs, Markdown reporting and signature/checksum receipts.
- Shared-server Conda environment, reference build script, schemas, examples and installation documentation.
- Existing STAR-index structural validation, provenance warnings and read-only coexistence auditing to avoid unnecessary index rebuilding.

## Validated in this workspace

- All Python source and tests compile under CPython 3.11.
- The dependency-free smoke test verifies REV plus/minus coordinate logic, reverse complementation, balanced-design acceptance, confounded-design rejection and deterministic contrast direction.
- Both JSON schemas parse successfully.
- Source contains no personal Windows/Linux home paths.
- R source has balanced delimiters; the Windows workspace has no R or Linux bioinformatics executables, so package-level R, Bash and end-to-end tests were not executable here.

## Required before production analysis

1. Build and lock the Linux Conda environment on the target server.
2. Populate and validate assembly-consistent human and/or mouse reference manifests.
3. Run the synthetic orientation/count/internal-priming truth set and at least one public QuantSeq REV benchmark.
4. Calibrate the initial PAS support, clustering and internal-priming thresholds rather than treating them as universal constants.
5. Pin PolyAseqTrap source/model assets and run its REV acceptance pilot. APA-B must remain disabled until that gate passes.
6. Review representative PCPA loci and use orthogonal assays before claiming actual premature transcription termination.

The workflow starts from demultiplexed `.fastq.gz` files. Instrument BCL/CBCL conversion and demultiplexing remain an upstream, facility-specific responsibility.
