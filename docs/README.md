# Documentation index

This index separates current user and administrator guidance from historical development records. For an ordinary analysis, start with the quick start and workflow-steps pages; historical alpha plans are not run instructions.

## Start and operate a project

| Question | Current document |
|---|---|
| How do I start an analysis on the shared server? | [Shared-server quick start](01_quick_start.md) |
| What does every workflow stage do? | [Workflow steps and dependencies](02_workflow_steps.md) |
| Which `config.conf` settings should I change? | [Configuration guide](03_configuration_guide.md) |
| How do I describe samples, replicates and paired libraries? | [Samplesheet contract](SAMPLESHEET_CONTRACT.md) |
| How do I monitor, resume or diagnose a run? | [Recovery and troubleshooting](recovery_and_troubleshooting.md) |
| How is the shared release installed? | [Server installation](server_installation.md) |

## Methods, QC and outputs

| Topic | Current document |
|---|---|
| Count universes and statistical methods | [Methods](methods.md) and [C0-C5 data stages](data_stages.md) |
| FastQ Screen species/contamination QC | [FastQ Screen](fastq_screen.md) |
| RSeQC orientation, feature distribution and gene-body coverage | [RSeQC](rseqc.md) |
| Track families, normalization and UCSC descriptors | [Tracks and outputs](tracks_and_outputs.md) |
| DGE/APA enrichment, plots and final reporting | [Enrichment and reporting](enrichment_and_reporting.md) |
| APA-B and comparison with APA-A | [APA-B and comparison](11_apa_b_and_comparison.md) |
| PAS atlas construction and provenance | [PAS atlases](pas_atlases.md) and [PAS atlas v1 design](PAS_ATLAS_V1_DESIGN.md) |
| Scientific and implementation limits | [Limitations](limitations.md) |

## Administrator and method-development references

- [PolyAseqTrap/DeepIP adapter and validation contract](POLYASEQTRAP_ADAPTER_CONTRACT.md)
- [PAS atlas v1 design](PAS_ATLAS_V1_DESIGN.md)

## Historical records

The following files explain earlier design and deployment decisions. They are retained for provenance but must not be used to install or run the current release:

- [Alpha.3 deployment-readiness review](ALPHA3_DEPLOYMENT_READINESS.md)
- [Alpha.1 GitHub/server deployment plan](GITHUB_AND_SERVER_DEPLOYMENT_PLAN.md)
- [Alpha.6 implementation snapshot](IMPLEMENTATION_STATUS.md)
- [Historical shared-server procedure](SHARED_SERVER_INSTALL.md)
- [Historical coexistence plan during CUT&RUN processing](SAFE_WORK_DURING_ACTIVE_CUTNRUN.md)

The current software version is recorded in the repository `VERSION` file and by `rna-ends2tracks --version`. Release changes belong in the root `CHANGELOG.md`.
