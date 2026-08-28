# Shared server installation

Install releases side-by-side into a writable shared environment parent. Do not modify an environment used by an active analysis.

```bash
git clone https://github.com/MichalGd/rna_ends2tracks.git
cd rna_ends2tracks
bash scripts/bash/install_release.sh \
  --tag v0.1.0-alpha.10.post6 \
  --env-parent /opt/conda_envs \
  --bin-dir /opt/conda_envs/bin \
  --mamba /opt/miniconda/condabin/mamba
```

The installer checks out the exact tag, creates a new environment, installs the package, runs Python tests, parses all R modules, executes the C4 paired/unpaired DESeq2 and track-factor smoke test, exports explicit/pip inventories, freezes the versioned environment, creates a self-contained versioned launcher and atomically promotes the stable symlink. The launcher prepends its immutable environment to `PATH` and clears inherited Python/R library overrides, so users do not need to activate Conda. Existing versioned launchers remain rollback targets.

Installation can run while another workflow is active because it creates different paths. Keep Mamba solving/installation at a low nice priority and choose conservative global run resources. Promotion affects only future invocations: an already running process continues with its loaded interpreter and absolute executable paths.

Reference FASTA/GTF/STAR/atlas assets remain shared read-only inputs and are not duplicated into the Conda environment. Validate both human and mouse config-only canaries before production, then run a small real-read preprocessing/exact-end canary for each genome.

## Optional APA-B environment

Install the pinned PolyAseqTrap/DeepIP branch separately after the main alpha.10 release:

```bash
bash scripts/bash/install_apa_b.sh --tag v0.1.0-alpha.10.post4
```

This creates a release-specific environment such as `/opt/conda_envs/rna_ends2tracks-apa-b-0.1.0a10.post4` without changing the stable main launcher. It is safe to install beside an active main-workflow run if disk and memory headroom are adequate. APA-B remains disabled until the synthetic pilot and a real QuantSeq REV canary for each intended assembly pass. Follow [the APA-B pilot contract](POLYASEQTRAP_ADAPTER_CONTRACT.md); do not set `APA_B_PILOT_ACCEPTED=true` merely because installation succeeded.

Rollback changes only the stable symlink:

```bash
ln -sfn /opt/conda_envs/bin/rna-ends2tracks-0.1.0-alpha.5 \
  /opt/conda_envs/bin/.rna-ends2tracks.rollback
mv -Tf /opt/conda_envs/bin/.rna-ends2tracks.rollback /opt/conda_envs/bin/rna-ends2tracks
```

Do not delete older environments until no project receipt or user command references them.
