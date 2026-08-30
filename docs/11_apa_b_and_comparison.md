# 11. APA-B and comparison with APA-A

APA-B answers the same biological question with a deliberately different discovery method. It clusters genome-wide QuantSeq REV 3′ endpoints with PolyAseqTrap, filters sequence-based internal-priming candidates with DeepIP, and tests PAS usage using DRIMSeq/stageR. APA-A remains the transparent Mcell2019-style primary analysis.

Both methods can detect PAS inside gene bodies. An internal-exonic or intronic PAS upstream of the normal terminal PAS is reported as a candidate premature cleavage/polyadenylation (PCPA) event. The assay does not directly measure polymerase occupancy, so the report must not call it proven premature transcription termination without independent evidence.

APA-B is implemented and available when the selected assembly and library protocol are covered by an accepted site validation manifest. The audited GRCm39 QuantSeq REV single-end new-project template enables both APA-A and APA-B by default. This default does not weaken the gate: an unvalidated assembly or protocol fails configuration validation unless APA-B is explicitly disabled or a matching accepted manifest is supplied.

The audited `biolserv` GRCm39 QuantSeq REV single-end deployment passed its synthetic and real-data pilots and may be enabled with the accepted site manifests documented in the [shared-server quick start](01_quick_start.md). GRCh38 and paired-end APA-B remain outside that accepted scope until their own real-data canaries and manifests pass. Site validation is an administrator operation, not something every project user repeats. After acceptance, the normal project still has only two inputs: `config.conf` and `samplesheet.csv`; users do not activate the APA-B environment or export its paths.

The two PAS catalogs are never combined. `08_apa_comparison/` matches nearby sites only after independent statistics and reports:

- sites detected by both methods;
- APA-A-only and APA-B-only sites;
- concordant and discordant directions of PAS-usage change;
- proximity between method-specific coordinates.

Agreement is supportive evidence, not a replacement for inspecting read coverage, replicate consistency, gene structure, and internal-priming evidence. See [the installation and pilot contract](POLYASEQTRAP_ADAPTER_CONTRACT.md) for exact pins, commands, outputs, and limitations.

## Untestable DRIMSeq/stageR hypotheses

DRIMSeq can legitimately return `NA` p-values for genes or PAS whose usage hypothesis cannot be estimated from the available counts. APA-B does not convert these values to zero, one, or statistical significance. Before stageR adjustment, it retains only genes with a finite screening p-value and at least two PAS with finite confirmation p-values. Excluded hypotheses remain `NA` in the complete result table. The supported stageR `allowNA=TRUE` policy is also enabled defensively.

Each contrast writes `<contrast>.na_audit.tsv`. This records the total and `NA` screening and confirmation tests, the numbers of genes and sites admitted to stageR, and the number of adjusted `NA` values. It also records the fixed policy `untestable hypotheses remain NA and cannot be significant`. A contrast with no testable stageR hypotheses fails with an actionable error instead of producing an empty success result.

## Paired-model numerical safeguards

Paired contrasts retain the resolved `~ subject + condition` design and use DRIMSeq's regression path rather than a one-way shortcut. Every contrast uses a deterministic seed. If, and only if, a multifactor fit fails with a recognized DRIMSeq numerical zero-pattern error, the complete fit is retried with DRIMSeq's documented `add_uniform=TRUE` option. This adds a reproducible sub-count perturbation to zeros for numerical fitting; it does not change the integer PAS matrix on disk, remove subjects, or replace the paired design with an unpaired model.

Each contrast writes `<contrast>.fit_audit.tsv`. `fit_policy=standard` means no retry was required. `fit_policy=deterministic_add_uniform_retry` records the fallback, seed, original error, and `WARN_NUMERIC_RETRY` status so the event remains visible in provenance and review.
