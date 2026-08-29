# Changelog

## Unreleased

## 0.1.0-alpha.10.post7 - 2026-08-29

- Overlap the DGE-then-final-tracks, APA-A, and APA-B branches under an aggregate thread/memory preflight ceiling; add separate APA-A and APA-B contrast-pool controls and publish final tracks without waiting for APA-B.
- Reconcile live status with successful stage receipts and count enrichment jobs from the authoritative enrichment index.
- Make cleanup safe across patch-release resumes and prior receipt-validated cleanup removals without weakening output-integrity checks.
- Automatically regenerate the exact APA-A statistical layer when a legacy index lacks gene-level enrichment sources.
- Expand the final scientific report with dedicated DGE and APA summaries, top differential genes, top APA gene events, method agreement percentages, and top enrichment terms.

## 0.1.0-alpha.10.post6 - 2026-08-28

- Keep paired APA-B contrasts on the resolved `~ subject + condition` model and explicitly use DRIMSeq's regression path for multifactor designs.
- Make DRIMSeq fitting reproducible with contrast-specific deterministic seeds.
- Retry only recognized multifactor numerical zero-pattern failures with DRIMSeq's documented `add_uniform=TRUE` safeguard, recording the original failure and selected policy in a per-contrast fit audit.

## 0.1.0-alpha.10.post5 - 2026-08-28

- Allow subset NA p-values from DRIMSeq to remain explicitly untestable during stageR adjustment instead of aborting every APA-B contrast.
- Exclude genes without a finite screening p-value and at least two finite PAS confirmation tests from stageR while preserving their output rows as `NA` and never significant.
- Record and validate a per-contrast APA-B NA audit covering screening, confirmation, stageR input, and adjusted-result counts.

## 0.1.0-alpha.10.post4 - 2026-08-28

- Reconstruct APA-B's unchanged raw endpoint universe from receipt-validated C1+C1S tables by default, avoiding a second full scan of every BAM while never reusing APA-A-specific C2 filtering.
- Retain an explicit BAM fallback and enforce C0=C1+C1S, source-receipt checksum, duplicate-retention, soft-clip, and endpoint-count invariants.
- Split endpoint preparation, PolyAseqTrap clustering, and DeepIP concurrency controls; use process-based endpoint preparation and bounded external R workers.
- Add resumable per-sample endpoint/clustering checkpoints plus live completed/total and ETA messages to the main engine log.
- Bound TensorFlow/BLAS CPU threads explicitly and record each sample's endpoint source in engine provenance.
- Pin the workflow-adapter commit in installation/acceptance provenance and require a new C1+C1S equivalence pilot.

## 0.1.0-alpha.10.post3 - 2026-08-28

- Resolve `Rscript` from the same immutable APA-B environment as the running adapter instead of relying on the caller's `PATH`.
- Remove inherited R library overrides from PolyAseqTrap subprocesses so both the synthetic pilot and production adapter use the pinned, installation-tested R package library.

## 0.1.0-alpha.10.post2 - 2026-08-28

- Add a reproducible synthetic APA-B coordinate, strand, count-conservation, duplicate-retention, intragenic-site, PolyAseqTrap-clustering, and pinned DeepIP truth-set pilot command.
- Add an explicit pre-acceptance adapter pilot mode for real QuantSeq canaries while preserving mandatory accepted-manifest validation for production execution.
- Install APA-B environments side-by-side under release-specific prefixes so immutable pilot environments can be upgraded without modification or deletion.
- Preserve pinned PolyAseqTrap `simpleCluster` ranges and centers while recalculating the intended per-cluster weighted scores from each cluster's `revmap`, preventing the upstream scalar global score from being recycled to every PAC.

## 0.1.0-alpha.10.post1 - 2026-08-28

- Resolve the APA-B environment's R/TensorFlow ICU conflict by solving the R 4.3 and Bioconductor stack with Conda while installing the pinned TensorFlow 2.10.1 CPU and Keras 2.10.0 wheels through pip.
- Verify the CPU-only DeepIP runtime before downloading and installing the pinned PolyAseqTrap and DeepIP sources.

## 0.1.0-alpha.10 - 2026-08-28

- Add embedded C4 PCA and sample-distance figures plus per-contrast MA and volcano plots in PDF/PNG.
- Add bounded, receipt-backed DGE/APA-A/PCPA and validated APA-B ORA/GSEA using GO, Reactome, and Hallmark gene sets with mapping/database provenance.
- Add a complete provenance dashboard covering configuration and samplesheet checksums, references/PAS atlases, stage receipts, environment packages, executable versions, and all output files.
- Require an accepted APA-B validation manifest and matching run-specific engine provenance before producing or interpreting APA-B results.
- Add a pinned, separately installable PolyAseqTrap/DeepIP APA-B adapter for QuantSeq REV: genome-wide weighted endpoint/PAC clustering, species-specific DeepIP filtering, intragenic PCPA eligibility, model/environment provenance, automated `config.conf` resolution, and synthetic/real-canary acceptance manifests.

- Run CPU-bound per-sample exact-end extraction and mixed Python/native track generation in bounded process pools so configured parallel jobs can use multiple CPU cores.
- Add an internal chronological master log, atomic run-status snapshot, and `rna-ends2tracks status` command.
- Record native-command and stage lifecycle events while retaining detailed per-tool logs.
- Explain the workflow-specific C0-C5 count universes in a beginner-facing data-stage guide.
- Split early C0 raw/CPM tracks from end-derived tracks with independent receipts and one strand-BAM extraction per sample/strand.
- Report completed/total work, elapsed time, and approximate ETA for bounded lane, sample, track, and contrast pools in the chronological master log.
- Start raw/CPM C0 tracks as each final sample BAM is published during merging, using only CPU/RAM capacity left after reserving the configured merge pool; defer automatically when no safe overlap capacity remains.
- Add a real tabular scientific HTML report and machine-readable per-contrast summary covering DGE direction, APA-A sites/shifts/PCPA, optional APA-B, and method concordance; abort reporting if source indexes disagree with result tables.
- Generate a one-column BigWig collection inventory, per-collection and combined one-line UCSC custom-track descriptors, distinct collection colors, configurable public URLs/minus-track display, and an IGV session during every report stage.
- Expand the HTML report with validated samples, STAR mapping, QuantSeq REV orientation, C0-C2R filtering funnels, active-PAS assignment, a searchable contrast table, track collection counts, and direct links to major QC and result indexes.

## 0.1.0-alpha.9.post2 - 2026-08-27

- Fix APA-A metadata mapping and fallback-comparator selection when DEXSeq returns factor-valued feature IDs: coerce PAS IDs to character and resolve them explicitly against the active-PAS catalog and named project C3 counts.
- Restore project PAS identifiers on the DEXSeq normalized-count matrix before computing mean counts and PAU, instead of matching against DEXSeq's internal composite row names.
- Make fallback-comparator filtering explicitly NA-safe so incomplete numerical rows cannot materialize artificial all-NA candidates during R data-frame subsetting.
- Compute pooled raw counts once per contrast and fail with an actionable message if a comparator PAS is absent from the C3 matrix.
- Add regression coverage for factor-valued PAS IDs, deterministic requested-ID order, missing PAS detection, and alpha.9.post1 receipt reuse.

## 0.1.0-alpha.9.post1 - 2026-08-27

- Fix unpaired APA-A DEXSeq formulas so a design without pairing covariates does not create a stray `:exon` term.
- Remove non-tabular list-valued DEXSeq metadata columns before writing result TSVs.
- Add a base-R paired/unpaired formula and result-serialization regression test that runs in CI and during installation.
- Allow immutable post-release tags so urgent alpha hotfixes do not consume the broader alpha.10 development version.
- Accept alpha.9 receipts only from this alpha.9.post1 hotfix, allowing downstream resume and final cleanup without recomputing successful upstream stages.

## 0.1.0-alpha.9 - 2026-08-27

- Run MultiQC without its incompatible cleanup inside a workflow-owned temporary directory, then
  remove that directory safely when package templates originate from an immutable shared environment.

## 0.1.0-alpha.8 - 2026-08-26

- Define the conservative GRCh38/v42 and GRCm39/vM31 PAS atlas v1 evidence policy.
- Add deterministic normalization of GENCODE polyA GTF and PolyA_DB v4.1 Main/Max archives.
- Require unique, strand-preserving mm10-to-GRCm39 conversion for mouse PolyA_DB sources.
- Add a fixed-source download helper with explicit UCSC chain-license acceptance.
- Replace the versioned launcher symlink with an immutable wrapper that exposes the release's own command-line and R tools without Conda activation.
- Clear inherited Python and R library overrides in the launcher and add a clean-environment regression test.

## 0.1.0-alpha.7 - 2026-08-26

- Fix the shared installer immutability audit so normal Conda symbolic links do not cause a false failure.
- Add a regression test that distinguishes writable symlinks from writable real files and directories.

## 0.1.0-alpha.6 - 2026-08-26

- Replaced the normal YAML interface with one restricted, commented `config.conf` plus one combined lane-level samplesheet.
- Added independent bounded pools for QC/trimming, STAR, merging, end extraction, tracks and pairwise statistical jobs, with global CPU/RAM preflight ceilings and timing records.
- Added simultaneous GRCh38/GRCm39 metadata support while strictly separating reference use, PAS discovery and contrasts by genome.
- Implemented C0-C5 count universes, end-defining clip auditing, refined internal-priming masking/rescue, and condition-blind two-round Mcell2019 PAS discovery.
- Made C4 active-PAS gene sums primary DGE and retained C5 featureCounts as an explicit diagnostic with correlation/discrepancy tables.
- Implemented deterministic DEXSeq comparator/shift rules, candidate intragenic PCPA reporting and pilot-gated independent APA-B comparison.
- Added the five default strand-specific track families, CPM for every family, and DESeq2/robust-CPM scaling for C2/C3.
- Added run locking, per-unit receipts, atomic outputs, success-only allow-list cleanup, PAS-atlas builder, side-by-side release installer and expanded documentation with Mermaid architecture.

## 0.1.0-alpha.5 - 2026-08-25

- **Breaking samplesheet contract:** existing sheets must add `description`, `genome`, and `technical_replicate_id`; active alpha.4 projects must finish with their versioned alpha.4 launcher.
- Added explicit `genome`, `technical_replicate_id`, and free-text `description` samplesheet columns.
- Normalize `hg38`/`GRCh38` and `mm39`/`GRCm39` aliases, require one genome per project, and reject a samplesheet/reference-manifest assembly mismatch before compute.
- Model biological samples, technical library preparations, and sequencing lanes separately. Multiple technical preparations and lanes sharing one `sample_id` are aligned independently and merged before DGE or APA statistics.
- Record technical-replicate and lane counts in validated biological-sample metadata and preserve technical-library identity in BAM read groups, filenames, logs, and orientation QC.
- Document the distinction between a biological replicate, technical replicate, sequencing lane, and matched-pair `subject`.

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
