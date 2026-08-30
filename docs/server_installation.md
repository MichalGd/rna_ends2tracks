# Shared server installation

Install releases side-by-side into a writable shared environment parent. Do not modify an environment used by an active analysis.

```bash
git clone https://github.com/MichalGd/rna_ends2tracks.git
cd rna_ends2tracks
bash scripts/bash/install_release.sh \
  --tag v0.1.0-alpha.11.post3 \
  --env-parent /opt/conda_envs \
  --bin-dir /opt/conda_envs/bin \
  --mamba /opt/miniconda/condabin/mamba
```

The installer checks out the exact tag, creates a new environment, installs the package, runs Python tests, parses all R modules, executes the C4 paired/unpaired DESeq2 and track-factor smoke test, exports explicit/pip inventories, freezes the versioned environment, creates a self-contained versioned launcher and atomically promotes the stable symlink. The launcher prepends its immutable environment to `PATH` and clears inherited Python/R library overrides, so users do not need to activate Conda. Existing versioned launchers remain rollback targets.

FastQ Screen itself and Bowtie2 are installed in the immutable workflow environment. Its large multi-genome indexes and `fastq_screen.conf` remain shared site resources outside that environment. Audit their readability, index/aligner compatibility, and reference provenance, then place the configuration path in each project's `FASTQ_SCREEN_CONFIG`. See [FastQ Screen QC](fastq_screen.md).

Installation can run while another workflow is active because it creates different paths. Keep Mamba solving/installation at a low nice priority and choose conservative global run resources. Promotion affects only future invocations: an already running process continues with its loaded interpreter and absolute executable paths.

Reference FASTA/GTF/STAR/atlas assets remain shared read-only inputs and are not duplicated into the Conda environment. Validate both human and mouse config-only canaries before production, then run a small real-read preprocessing/exact-end canary for each genome.

## Separately validated APA-B environment

Install the pinned PolyAseqTrap/DeepIP adapter in its separate immutable environment:

```bash
bash scripts/bash/install_apa_b.sh --tag v0.1.0-alpha.10.post6
```

This creates `/opt/conda_envs/rna_ends2tracks-apa-b-0.1.0a10.post6` without changing the stable main launcher. The current main workflow resolves this sidecar through its installation manifest; users never activate it manually. It is safe to install beside an active main-workflow run if disk and memory headroom are adequate.

Installation and scientific acceptance are separate. The audited `biolserv` GRCh38/GRCm39 QuantSeq REV V2 single-end deployment has passed the synthetic and assembly-specific real-data pilots and may use its combined accepted validation manifest. Each additional assembly or library protocol—particularly paired-end APA-B—requires a matching real-data canary and accepted manifest. Follow [the APA-B pilot contract](POLYASEQTRAP_ADAPTER_CONTRACT.md); do not set `APA_B_PILOT_ACCEPTED=true` merely because installation succeeded.

Rollback changes only the stable symlink:

```bash
ln -sfn /opt/conda_envs/bin/rna-ends2tracks-0.1.0-alpha.5 \
  /opt/conda_envs/bin/.rna-ends2tracks.rollback
mv -Tf /opt/conda_envs/bin/.rna-ends2tracks.rollback /opt/conda_envs/bin/rna-ends2tracks
```

Do not delete older environments until no project receipt or user command references them.
