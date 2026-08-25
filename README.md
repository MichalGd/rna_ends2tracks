# rna_ends2tracks

`rna_ends2tracks` is a configuration-driven command-line workflow for Lexogen QuantSeq REV libraries without UMIs. It performs shared preprocessing/alignment, independent gene-expression analysis, two analytically independent APA analyses, intragenic premature-cleavage classification, browser-track generation and cross-method reporting. It deliberately uses standalone Python, Bash and R modules rather than a workflow engine.

The same installation analyses human GRCh38 and mouse GRCm39. Species-specific behavior is confined to validated reference manifests; no user-specific or server-specific path is embedded in the code.

## Implemented analysis branches

- Shared per-lane FastQC, Lexogen-aware BBDuk adapter/poly(A)/poly(T)/quality trimming, post-trim QC, lane-specific STAR read groups, technical-lane merging, samtools validation, MultiQC and empirical strandedness checks.
- Independent exon-union featureCounts (`-s 2`) and DESeq2 analysis for every eligible condition pair.
- APA-A exact QuantSeq REV R1 end extraction, unique-primary filtering, soft-clip audit, condition-blind clustering, site-level internal-priming evidence, GTF annotation, raw counts, DEXSeq/Delta-PAU and statistically filtered PCPA candidates.
- Pilot-gated APA-B adapter for an independently pinned PolyAseqTrap installation, output-contract validation, DRIMSeq/stageR, and independently classified PCPA candidates.
- Strand-aware APA crosswalk, effect-direction comparison, PCPA agreement and transcript-strand BigWigs.

Coordinate-identical reads and BAM duplicate flags are retained throughout. UDI sample indexes are never treated as UMIs.

## Requirements and installation

Use Linux on a shared server or HPC system. Create the supplied Conda environment and install the CLI:

```bash
mamba env create -f environment.yml
conda activate rna_ends2tracks-0.1.0a3
python -m pip install --no-deps .
rna-ends2tracks --version
```

The source installation above is appropriate for development or an alpha canary. A shared release should install the tagged wheel with `python -m pip install --no-deps <wheel>`; do not use `pip install -e` in a shared production environment.

See [shared-server installation](docs/SHARED_SERVER_INSTALL.md) for a read-only installation usable by every server user. PolyAseqTrap is isolated behind the [APA-B adapter contract](docs/POLYASEQTRAP_ADAPTER_CONTRACT.md) because its REV behavior and DeepIP assets must pass a local pilot before use.

Current verification and production-gate details are recorded in [implementation status](docs/IMPLEMENTATION_STATUS.md).

Repository publication, side-by-side installation, concurrent-analysis safety and rollback are covered by the [GitHub and shared-server deployment plan](docs/GITHUB_AND_SERVER_DEPLOYMENT_PLAN.md).

For deployment while `cutnrun2tracks` jobs are active, follow the stricter [active-CUT coexistence runbook](docs/SAFE_WORK_DURING_ACTIVE_CUTNRUN.md) and begin with its read-only server audit.

## Configure a project

1. Copy `config/project.example.yaml` and `config/samplesheet.example.csv` into the project.
2. Copy either `references/human_GRCh38.example.yaml` or `references/mouse_GRCm39.example.yaml` and replace every placeholder with a pinned, assembly-consistent asset.
3. Set `condition_order`. For each unordered pair of conditions having at least two biological replicates, the later condition is the numerator and the earlier condition the denominator.
4. Use one samplesheet row per sequencing lane. Repeated rows for a `sample_id` are technical lanes, never statistical replicates.
5. Keep `umi_present: false`, `protocol.has_umi: false`, and `retain_duplicate_flagged_reads: true`.

Only validated single-end QuantSeq REV profiles are enabled. Paired-end REV is rejected until its R2 first-12-base behavior and end-coordinate contract receive a separate pilot. The design formula currently supports additive terms such as `~ batch + condition` or `~ subject + condition`; rank-deficient/confounded designs fail before compute.

## Run

```bash
rna-ends2tracks --config project.yaml --samplesheet samples.csv validate
rna-ends2tracks --config project.yaml --samplesheet samples.csv --dry-run all
rna-ends2tracks --config project.yaml --samplesheet samples.csv all
```

Modules can be scheduled independently:

```bash
rna-ends2tracks --config project.yaml --samplesheet samples.csv preprocess
rna-ends2tracks --config project.yaml --samplesheet samples.csv dge
rna-ends2tracks --config project.yaml --samplesheet samples.csv apa-a
rna-ends2tracks --config project.yaml --samplesheet samples.csv tracks
```

APA-B remains off in the example. After the documented pilot succeeds, configure its adapter, enable it and run:

```bash
rna-ends2tracks --config project.yaml --samplesheet samples.csv apa-b
rna-ends2tracks --config project.yaml --samplesheet samples.csv compare
rna-ends2tracks --config project.yaml --samplesheet samples.csv report
```

The source-tree launcher `bash bin/run_module.sh ...` is equivalent when the package is not installed.

## Output structure

```text
results/
├── 00_metadata/       validated samples, lanes, contrasts, resolved configuration
├── 01_qc/             raw/trimmed FastQC, MultiQC and trimmed FASTQs
├── 02_alignment/      lane and biological-sample BAMs, indexes and protocol audit
├── 03_gene_expression/
├── 04_apa_a_repository/
├── 05_apa_b_polyaseqtrap_drimseq/
├── 06_apa_comparison/
├── 07_tracks/
├── 08_reports/
└── provenance/        commands, structured events and completion receipts
```

A module is considered complete only when its validated `run_receipt.json` exists and its signature and output checksums still match. `--force-module` reruns one selected module; the workflow has no recursive cleanup command.

## Scientific limitation

An intragenic QuantSeq end can support a candidate premature cleavage and polyadenylation event. It cannot by itself prove premature transcription termination. Reports therefore use “candidate PCPA consistent with premature transcription termination.” Claims of downstream polymerase loss require nascent-transcription evidence; transcript structure should be confirmed with long reads, 3′ RACE or targeted assays.
