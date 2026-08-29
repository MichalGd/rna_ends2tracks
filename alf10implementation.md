# Alpha.10 APA-B implementation record

Date: 2026-08-28  
Development branch: `feature/alpha10-development`  
Primary implementation commit: `05165f3` (`Implement pinned QuantSeq APA-B adapter`)

## Outcome

Alpha.10 implements an independent second APA analysis path, APA-B, based on PolyAseqTrap clustering and DeepIP internal-priming classification. It is designed for non-UMI Lexogen QuantSeq REV data and can retain cleavage sites throughout annotated genes, including intronic and internal-exonic candidate premature cleavage and polyadenylation sites far upstream of annotated transcript ends.

The implementation includes:

- genome-wide, strand-aware transcript-end extraction from C0 BAM files;
- weighted PolyAseqTrap `simpleCluster` clustering of cleavage endpoints;
- project-wide, condition-blind merging and support filtering of candidate sites;
- species-specific DeepIP classification for human and mouse;
- an APA-B count matrix and independent DRIMSeq/stageR differential-usage analysis;
- annotation of candidate intragenic premature cleavage and polyadenylation events;
- an immutable, separately installed APA-B environment with pinned sources and checksums;
- synthetic and real-data pilot validation machinery;
- post-hoc APA-A/APA-B proximity and effect-concordance comparison without merging their catalogs.

## Post6 production-readiness follow-up

The first complete 18-sample run exposed operational issues outside the APA-B scientific model. Release post7 adds budget-checked concurrent execution of the DGE-then-final-tracks, APA-A, and APA-B branches; separate APA-method contrast pools; correct enrichment/status accounting; safe cross-patch cleanup and repeat-cleanup evidence handling; exact legacy APA-A statistical regeneration before enrichment; and expanded DGE, APA, concordance, top-event, and enrichment summaries in the final report. Final tracks can therefore publish while APA-B is still running. These changes do not merge APA catalogs or alter either method's statistical universe.

## QuantSeq-specific scientific decision

Unmodified PolyAseqTrap `FindPTA` prioritizes reads with primary poly(A)-tail evidence and stops when that evidence is absent. Lexogen QuantSeq REV libraries do not reliably retain a sequenced poly(A) tail, and the main workflow trims terminal poly(A/T) sequence before STAR alignment.

For that reason, APA-B is an explicitly documented assay adaptation:

1. extract strand-aware genomic transcript 3′ endpoints from aligned QuantSeq reads;
2. cluster weighted endpoints using the exported PolyAseqTrap `simpleCluster` algorithm with a 24-nt clustering distance;
3. form a project-wide, condition-blind PAS universe;
4. apply reproducible sample/read-support filters;
5. classify candidate sites with the official species-specific DeepIP model;
6. quantify retained sites and test differential usage independently of APA-A.

The workflow records this as `polyaseqtrap_quantseq_rev`, not as unmodified `FindPTA` execution.

## Independence from APA-A

APA-B shares only necessary upstream biological inputs with APA-A: C0 BAM files, the reference FASTA, and the reference annotation. It does not consume APA-A active-PAS catalogs, PAS-atlas rescue decisions, site identifiers, internal-priming calls, or APA-A statistical thresholds.

APA-A and APA-B outputs remain separate. Their comparison is a downstream proximity/effect-concordance analysis and cannot alter either source catalog.

## Main implementation files

- `src/rnaends2tracks/polyaseqtrap_adapter.py` — endpoint extraction, orchestration, validation, and output contract.
- `scripts/R/polyaseqtrap_quantseq_rev.R` — PolyAseqTrap clustering, DeepIP integration, annotation, and count generation.
- `environment.apa_b.yml` — separate APA-B software environment.
- `scripts/bash/install_apa_b.sh` — pinned immutable installation and manifest creation.
- `src/rnaends2tracks/apa_b_pilot.py` — validation evidence and accepted-manifest builder.
- `docs/POLYASEQTRAP_ADAPTER_CONTRACT.md` — exact adapter contract and scientific limitations.
- `docs/11_apa_b_and_comparison.md` — user-facing APA-B operation and interpretation.
- `tests/test_polyaseqtrap_adapter.py` — adapter and contract regression tests.

## Pinned upstream components

- PolyAseqTrap commit: `176ea2884ff1c6be7c64bc44fa7661d82d90e718`
- DeepIP commit: `988564875d002b6d5d48d8dfb228cba3492dd776`
- Human DeepIP model SHA-256: `d74138c788102ae57a50664b6858a0b79951430fee9bdcc93b07f9b1ba16edf1`
- Mouse DeepIP model SHA-256: `aba432c85ef6c14e56a6222106acaffbcc3b9131a86508afdf66311fe57123e9`

Upstream references:

- [PolyAseqTrap](https://github.com/APAexplorer/PolyAseqTrap)
- [PolyAseqTrap clustering and FindPTA source](https://github.com/APAexplorer/PolyAseqTrap/blob/176ea2884ff1c6be7c64bc44fa7661d82d90e718/R/PolyAseqTrap_funclib.R#L99-L286)
- [DeepIP](https://github.com/APAexplorer/DeepIP)
- [DeepIP inference script](https://github.com/APAexplorer/DeepIP/blob/988564875d002b6d5d48d8dfb228cba3492dd776/DeepIP_test.py)

## Configuration boundary

APA-B should remain disabled until the separate APA-B installation and the site-specific pilot have passed. After an accepted validation manifest is created, the relevant configuration is:

```bash
RUN_APA_B=true
APA_B_PILOT_ACCEPTED=true
APA_B_COMMAND_TEMPLATE="auto"
APA_B_INSTALLATION_MANIFEST="/opt/conda_envs/rna_ends2tracks-apa-b-0.1.0a10.post4/installation_manifest.json"
APA_B_VALIDATION_MANIFEST="/opt/conda_envs/rna_ends2tracks-apa-b-validation/accepted_GRCm39.json"
APA_B_THREADS=8
```

This gate means that the implementation is present but production scientific use is deliberately blocked until the installed models and real QuantSeq behavior have been validated on the server.

## Validation completed in development

The standard GitHub CI passed for the implementation commit, including:

- Python 3.10 and 3.11 tests;
- package-asset and configuration-schema checks;
- Bash syntax and ShellCheck;
- R parsing and repository-contract tests.

At the time of this record, the complete-environment GitHub job and the server-side APA-B scientific pilot are separate remaining release gates. The already completed alpha.9 production analysis is not affected by this branch.

## Remaining release and server steps

1. Confirm every GitHub check on `feature/alpha10-development` is green and merge it into `main`.
2. Tag the merged commit as the alpha.10 release.
3. Install the separate APA-B environment on the server:

   ```bash
   bash scripts/bash/install_apa_b.sh --tag v0.1.0-alpha.10.post4
   ```

4. Run the synthetic pilot, verifying coordinate, strand, count, clustering, and DeepIP invariants.
5. Run small real human and mouse QuantSeq REV canaries.
6. Review the evidence and create the accepted validation manifest; the workflow must not self-approve it.
7. Enable APA-B in `config.conf` and resume from the APA-B step.

## Interpretation limits

- An intragenic APA-B site is a candidate premature cleavage and polyadenylation site; 3′ RNA-seq alone does not prove RNA polymerase II transcription termination.
- APA-A/APA-B agreement increases confidence but is not a substitute for orthogonal validation.
- DeepIP model applicability and classification behavior must be evaluated separately for each supported organism and reference build.
- PolyAseqTrap is GPL-3.0. The DeepIP repository did not expose an unambiguous redistribution license during implementation review, so institutional approval or author clarification is required before redistributing its code or model files.
