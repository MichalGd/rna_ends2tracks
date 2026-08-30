# PolyAseqTrap/DeepIP adapter and QuantSeq REV pilot contract

## What is implemented

APA-B is a repository-owned, installable alternative APA method. It uses pinned versions of:

- PolyAseqTrap commit `176ea2884ff1c6be7c64bc44fa7661d82d90e718`;
- DeepIP commit `988564875d002b6d5d48d8dfb228cba3492dd776`;
- the official human or mouse DeepIP model, verified by SHA-256;
- a separate explicit Conda environment lock.

The adapter is independent from APA-A: it never reads the Mcell active-PAS catalog, PAS-atlas rescue decisions, APA-A identifiers, or APA-A internal-priming calls. Only the common C0 alignment evidence, assembly FASTA, and GTF are shared. APA-A and APA-B are compared after both analyses finish.

## Why a QuantSeq-specific adaptation is necessary

PolyAseqTrap `FindPTA()` prioritizes reads that retain informative poly(A)-tail sequence and stops when no primary V1/V2 evidence is present. Lexogen QuantSeq REV reads do not reliably retain a sequenced poly(A) tail, and `rna_ends2tracks` removes technical poly(A/T) sequence before STAR. Calling unmodified `FindPTA()` on those BAMs would therefore be unreliable.

The adapter instead:

1. extracts every mapped primary C0 transcript 3′ endpoint without UMI or coordinate deduplication;
2. collapses identical endpoints to integer weights while retaining their full read count;
3. calls PolyAseqTrap `simpleCluster()` with its 24-nt weighted PAC clustering rule;
4. merges per-sample PACs project-wide, condition-blind, with deterministic strand-aware summit selection;
5. retains sites supported by at least five reads in at least two samples by default;
6. sends A-rich/repeat candidates to the official species-specific DeepIP model;
7. annotates terminal, internal-exonic, intronic, ambiguous, and intergenic sites against the assembly-matched GTF;
8. runs DRIMSeq/stageR and separately classifies candidate intragenic PCPA events.

The pinned `simpleCluster()` implementation supplies the official strand-aware cluster ranges and centers. Its scalar global score assignment is not used: the adapter recalculates each PAC's weighted count from that cluster's own `revmap` members and requires exact record-count conservation in the synthetic pilot.

Because discovery is genome-wide, an intronic or internal-exonic site far upstream of the annotated gene end is eligible. Such a QuantSeq signal is a **candidate premature cleavage/polyadenylation event**, not by itself proof of RNA polymerase II termination.

## Installation

APA-B has older deep-learning dependencies and is deliberately isolated from the main workflow environment:

```bash
bash scripts/bash/install_apa_b.sh --tag v0.1.0-alpha.10.post6
```

The default installation is release-specific, for example `/opt/conda_envs/rna_ends2tracks-apa-b-0.1.0a10.post6`. The script verifies the source/model pins, records the exact workflow-adapter release and Git commit, installs the adapter, exports an explicit environment lock, writes `installation_manifest.json`, and makes the versioned environment read-only. Installation alone does **not** accept the method scientifically.

PolyAseqTrap declares GPL-3. The reviewed DeepIP repository does not currently expose a clear software/model license file. Non-commercial research use does not itself resolve missing license terms; the server administrator should record institutional approval or obtain clarification from the authors before redistribution or broader deployment.

## Pilot gate

Before production use, run a synthetic truth set and at least one real QuantSeq REV canary for every assembly and protocol scope that will be enabled. The synthetic audit must confirm strand/base coordinates, exact equivalence of BAM-derived endpoints and receipt-validated C1+C1S reuse, count conservation, duplicate-flag retention, DeepIP positive/negative behavior, and retention of an intragenic site. A real canary must produce the five adapter deliverables and conserve eligible records. The audited `biolserv` post6 sidecar, synthetic pilot, and GRCh38/GRCm39 QuantSeq REV V2 single-end canaries satisfy this gate for those two assemblies and that protocol. Paired-end APA-B still needs its own accepted real-data validation.

Run the synthetic pilot directly from the immutable APA-B environment:

```bash
APA_B_ENV=/opt/conda_envs/rna_ends2tracks-apa-b-0.1.0a10.post6
"$APA_B_ENV/bin/rna-ends2tracks-run-apa-b-synthetic-pilot" --installation-manifest "$APA_B_ENV/installation_manifest.json" --output /path/to/synthetic_pilot_audit.json
```

For a pre-acceptance real canary, invoke `rna-ends2tracks-apa-b` with `--pilot-mode`. This mode verifies all installed engine/model checksums and records `pilot_mode=true` in provenance, but does not accept the method. Normal workflow execution never supplies this flag and still requires an accepted validation manifest.

Create the acceptance manifest only after reviewing those results:

```bash
/opt/conda_envs/rna_ends2tracks-apa-b-0.1.0a10.post6/bin/rna-ends2tracks-accept-apa-b-pilot \
  --installation-manifest /opt/conda_envs/rna_ends2tracks-apa-b-0.1.0a10.post6/installation_manifest.json \
  --synthetic-audit /path/to/synthetic_pilot_audit.json \
  --real-canary GRCh38=/path/to/human_canary/output \
  --real-canary GRCm39=/path/to/mouse_canary/output \
  --reviewed-by "REVIEWER NAME" \
  --output /opt/conda_envs/rna_ends2tracks-apa-b-validation/accepted_GRCh38_GRCm39.json
```

The builder refuses incomplete audits. Human and mouse may be accepted together by repeating `--real-canary`.

## Project configuration after acceptance

```bash
RUN_APA_B=true
APA_B_PILOT_ACCEPTED=true
APA_B_COMMAND_TEMPLATE="auto"
APA_B_INSTALLATION_MANIFEST="/opt/conda_envs/rna_ends2tracks-apa-b-0.1.0a10.post6/installation_manifest.json"
APA_B_ENDPOINT_SOURCE="auto"
APA_B_THREADS=16
APA_B_ENDPOINT_PARALLEL_JOBS=8
APA_B_CLUSTER_PARALLEL_JOBS=12
APA_B_DEEPIP_THREADS=16
APA_B_VALIDATION_MANIFEST="/opt/conda_envs/rna_ends2tracks-apa-b-validation/accepted_GRCh38_GRCm39.json"
```

`APA_B_ENDPOINT_SOURCE="auto"` first validates the existing per-sample exact-end receipt, checks the recorded counts and source files, and reconstructs the raw adapter input as C1+C1S. This is count- and coordinate-equivalent to rereading C0: C1 contributes unclipped exact ends and C1S contributes the separately reported end-clipped records that APA-B retains. C2 is never reused because it contains APA-A-specific internal-priming filtering. If a complete validated C1+C1S set is absent, `auto` falls back to BAM extraction; `exact_ends` makes reuse mandatory and `bam` forces an independent rescan.

Endpoint preparation uses processes, PolyAseqTrap uses bounded per-sample R subprocesses, and DeepIP receives explicit TensorFlow/BLAS thread limits. Successful endpoint and clustering work units receive checkpoints, so an interrupted adapter run can resume them. Progress and approximate ETA are written to the APA-B engine log.

An existing completed APA-A run can then resume at APA-B:

```bash
rna-ends2tracks --from-step apa_b /path/to/config/config.conf
```

## Required adapter outputs

For each assembly under `07_apa_b/`, the adapter creates:

- `pas_catalog.tsv` — independent genome-wide PAS/PAC coordinates, gene/feature annotation, support, and method evidence;
- `pas_counts.tsv` — integer raw PAC counts in exact validated sample order;
- `deepip_audit.tsv` — every evaluated candidate, prediction/score, and retention decision;
- `engine_provenance.json` — verified engine, model, environment, assay adaptation, and count-conservation evidence;
- `adapter_audit.json` — per-sample input/endpoint/duplicate/soft-clip accounting.

The workflow validates these files before differential analysis and refuses APA-A identifiers in the APA-B catalog or matrix.

## Scientific limitations

- The method is an assay-specific adaptation of PolyAseqTrap clustering and DeepIP, not an unmodified `FindPTA()` analysis.
- DeepIP resolves sequence-based internal-priming risk; it cannot prove that every retained site is a functional PAS.
- GTF assignment can be ambiguous where same-strand genes overlap; ambiguous sites remain visible but are excluded from gene-level statistics.
- Adding samples changes the condition-blind discovery universe. Use a new output directory.
- APA-A/APA-B agreement increases confidence but disagreement is not automatically an error; the methods intentionally use different discovery and filtering rules.
