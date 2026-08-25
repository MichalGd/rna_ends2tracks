# Safe work while two cutnrun2tracks analyses are running

## Current situation and principle

Two `cutnrun2tracks` instances are active. The public repository was inspected at commit `ea89ff2fbca4273adf0d47ac13451fc234470578` on 2026-08-25. Its documented environment is a useful package/reference baseline, but the public repository does not prove which commit, environment prefix or references the two server jobs currently use.

The safe rule is:

> Read and inventory the running installation, but publish `rna_ends2tracks` only into new versioned paths. Never modify anything resolved from an active process.

This permits useful work immediately without waiting for the CUT&RUN analyses to finish.

## Resources likely reusable after read-only validation

The documented CUT installation already includes many tools needed by `rna_ends2tracks`:

- Python 3.11;
- FastQC and MultiQC;
- samtools and bedtools;
- deepTools and `bedGraphToBigWig`;
- R and DESeq2;
- human GRCh38 and mouse GRCm39 FASTAs;
- GRCh38/GRCm39 chromosome-size files; and
- human/mouse GTF annotations.

Reuse has two meanings:

1. **References:** reuse validated FASTA, GTF and chromosome-size files directly through read-only absolute paths. This avoids copying large assets.
2. **Software baseline:** reuse known package versions, the local Conda package cache, or a `--copy` clone into a new environment. Do not make `rna_ends2tracks` call executables from the live CUT environment individually; cross-prefix binaries and libraries can be fragile.

The CUT Bowtie2 indexes cannot replace STAR indexes. A new STAR index is required unless a matching, validated STAR index already exists. A PAS atlas is also expected to be new unless one is found and assembly/contig compatibility is proven.

An existing STAR index is preferred over rebuilding when all of the following hold:

- `Genome`, `SA`, `SAindex`, `chrName.txt`, `chrLength.txt` and `genomeParameters.txt` are complete and readable;
- contig names and lengths exactly match the selected FASTA `.fai`;
- an existing manifest or checksums link the index to the same FASTA sequence and intended GTF release;
- it is an ordinary host index rather than an unintended host-plus-spike-in/composite index;
- the recorded STAR version is compatible with the runtime STAR;
- `sjdbOverhang` is suitable for the QuantSeq read length, or any mismatch has been scientifically reviewed;
- the annotation release is appropriate for gene counting and splice-junction discovery; and
- the path is immutable/read-only to workflow users.

Concurrent jobs may read the same immutable STAR index without corrupting it. The remaining concern is initial shared-filesystem I/O when several aligners load a large index at once. Schedule the RNA alignment separately if the two CUT jobs are currently I/O-heavy.

Reject or rebuild an index when contig names/lengths or recorded FASTA/GTF checksums differ, provenance is irrecoverable, required files are incomplete, the index is mutable, or it represents the wrong/composite assembly. An undocumented `sjdbOverhang` or different annotation is a review condition rather than an automatic rejection when all structural and sequence evidence agrees.

Likely missing from the CUT environment and required in the independent RNA environment:

- STAR;
- BBMap/BBDuk;
- Subread/featureCounts;
- PyYAML and pysam;
- DEXSeq;
- DRIMSeq and stageR; and
- later, the separately pinned PolyAseqTrap/DeepIP environment.

The provided audit script tests rather than assumes these observations.

## Actions already safe and implemented locally

- GitHub-ready repository boundary under `rna_ends2tracks/`.
- `.gitattributes`, version, changelog, citation metadata, contributing and security policies.
- GitHub Actions for Python tests/builds, installed package-data checks, Bash/ShellCheck, R parsing and repository hygiene.
- Explicit license-decision gate before public release.
- Read-only server coexistence audit script: `scripts/bash/audit_existing_server.sh`.

None of these actions contact or alter the running server jobs.

## Safe server work right now

### Step 1: snapshot the active jobs without signalling them

Use the site's scheduler first, for example `squeue -u "$USER"` on SLURM. A process snapshot is also read-only:

```bash
pgrep -af 'cutnrun2tracks|preprocess_batch|align_batch|peakcall_batch|differential_batch|coverage_batch'
```

Record each job ID, process ID, start time, current stage, output directory, workflow checkout, Conda prefix and reference paths. Read scheduler metadata, logs and `/proc/<pid>/cwd` only when permitted. Do not use `kill`, `renice`, debuggers, tracing or cleanup commands.

### Step 2: copy only the audit script to a user-owned staging location

Copy the script to a directory unrelated to both active output trees and the live CUT installation. Do not put it under the CUT `current` release or environment.

Run it with the observed environment and explicit known references:

```bash
bash audit_existing_server.sh \
  --cut-env /actual/read-only/cutnrun/environment \
  --output /user/project/rna_ends2tracks_audit_20260825 \
  --human-fasta /actual/path/GRCh38.fa \
  --human-gtf /actual/path/GRCh38.annotation.gtf \
  --human-chrom-sizes /actual/path/GRCh38.chrom.sizes \
  --human-star-index /actual/path/GRCh38_STAR_index \
  --human-pas-atlas /actual/path/GRCh38_PAS_atlas.bed.gz \
  --mouse-fasta /actual/path/GRCm39.fa \
  --mouse-gtf /actual/path/GRCm39.annotation.gtf \
  --mouse-chrom-sizes /actual/path/GRCm39.chrom.sizes \
  --mouse-star-index /actual/path/GRCm39_STAR_index \
  --mouse-pas-atlas /actual/path/GRCm39_PAS_atlas.bed.gz
```

The audit intentionally does not hash large references, load a STAR genome into memory, activate/clone/update an environment, create an index, modify permissions or signal jobs. It compares STAR contig names/lengths with an existing FASTA `.fai` and records index parameters. Its output is an inventory, not an authorization to reuse an asset.

### Step 3: inspect capacity

Read-only checks:

```bash
uptime
free -h
df -h /shared/apps /shared/references /project/scratch
```

Use scheduler/accounting commands preferred by the site. Avoid repeated broad filesystem scans. If the two CUT jobs are alignment, coverage, peak-calling or metagene stages, assume shared storage and CPU may already be busy.

### Step 4: prepare source and release artifacts off production paths

Safe options:

- push the reviewed private GitHub repository from the local workstation after owner/visibility/license decisions;
- build the Python wheel locally or in GitHub Actions;
- download the immutable source archive and wheel to a new user-owned staging directory;
- verify release SHA-256 checksums; and
- draft reference manifests using audited existing paths.

These steps do not require executing bioinformatics tools on the server.

### Step 5: optionally create the new environment, subject to load

Creating a new prefix does not change the active CUT environment, so it is correctness-safe. It can still contend for filesystem I/O and package-cache locks.

Proceed now only if:

- the new prefix is unused and outside every active environment/output directory;
- there is sufficient free disk space;
- the site permits environment construction on the chosen node;
- the work runs through a low-resource scheduler job or on a dedicated build node; and
- no command contains `conda update`, `pip install`, or `pip uninstall` against the live CUT prefix.

Two possible approaches:

**Preferred reproducible approach:** create a fresh prefix from a Linux lock file. This can reuse packages already present in the Conda package cache without inheriting CUT-specific state.

**Pilot reuse approach:** clone the confirmed CUT environment with `conda create --copy --clone` into a new RNA-specific prefix, record the explicit package set, then add only missing RNA packages with installed packages frozen. Abort if the solver proposes replacing the baseline R/Python/Bioconductor stack. This uses more I/O and should be scheduled.

Do not make the new environment visible as a default module yet.

### Step 6: run only lightweight canaries

While CUT analyses are active, acceptable checks on a low-load/build node are:

- `rna-ends2tracks --version`;
- Python unit tests;
- Bash syntax and R parsing;
- metadata-only human/mouse validation using `--skip-input-checks`; and
- checking package imports and executable availability.

Do not run STAR index generation, whole-reference hashing, real FASTQ processing, MultiQC scans over shared CUT outputs, or public-data benchmarks on the same busy server resources yet.

## Work that should be scheduled but need not wait for CUT completion

If the scheduler can isolate resources, submit these as separate low-priority jobs with new output paths:

- fresh Conda environment construction;
- wheel installation and complete unit tests;
- STAR index construction in a new versioned reference directory;
- tiny synthetic human and mouse integration tests; and
- one small QuantSeq pilot.

Specify CPU, RAM, walltime and scratch explicitly. The jobs must not depend on, scan, or write to either CUT project directory. Scheduler isolation prevents compute oversubscription, but shared filesystem I/O may still justify postponing STAR indexing until the CUT alignment/coverage stages finish.

## Absolutely defer

- Changing the CUT environment or its packages.
- Editing the checked-out CUT workflow used by either job.
- Replacing any reference or index named in an active command/configuration.
- Changing the live CUT or RNA `current` symlink.
- Making the alpha RNA module the default.
- Deleting Conda caches, old environments, old releases or references.
- Reusing either running job's scratch/output directory.
- Running automatic cleanup in either workflow.
- Enabling APA-B before the independent REV pilot.

## Go/no-go matrix for today

| Action | Correctness risk to CUT jobs | Resource interference | Decision now |
|---|---:|---:|---|
| Prepare/commit local Git repository | None | None | Go |
| Create private GitHub remote | None to server | Network only | Go after owner/visibility/license decision |
| Read scheduler/process/log metadata | None | Negligible | Go |
| Run provided read-only audit | None | Low | Go |
| Draft manifests pointing to existing references | None | None | Go |
| Build wheel in GitHub/local workstation | None | None | Go |
| Create fresh Conda prefix on dedicated build node | None if path is new | Moderate I/O | Conditional go |
| Clone CUT environment with `--copy` | None if source remains read-only | High metadata/I/O | Schedule/conditional |
| Build STAR index on isolated compute/scratch | None if output is new | High CPU/RAM/I/O | Scheduled conditional |
| Modify active CUT environment/reference | High | Variable | No-go |
| Publish/change `current` or default module | Can alter later job steps/new jobs | Low | No-go |
| Run real RNA data | None if fully isolated | High | Wait for canary and capacity |

## Minimum information still needed for server execution

- SSH hostname or deployment mechanism and an account authorized for the chosen staging/application paths.
- Scheduler type and current CUT job IDs/stages.
- Exact active CUT workflow path and environment prefix.
- Actual GRCh38/GRCm39 FASTA, GTF and chromosome-size paths.
- Existing STAR index and PAS atlas paths, if any.
- Site application/reference roots and permissions policy.
- Available disk/scratch and whether environment/index builds must be scheduled.
- GitHub repository owner, initial visibility, and explicit license choice.

Until these are known, local repository preparation and the read-only audit are the maximum safe implementation; guessing server paths would create the disturbance this plan is designed to avoid.
