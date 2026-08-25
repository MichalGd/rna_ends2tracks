# Shared-server installation

The workflow contains no usernames, home-directory assumptions or fixed project paths. Install one read-only software environment and reference set, while each user writes results to a project directory they own.

## Recommended administrator procedure

1. Clone or copy a tagged workflow release under a versioned, read-only path such as `/opt/rna_ends2tracks/0.1.0`.
2. Create the Conda environment from `environment.yml` with Mamba/Conda. Generate a platform-specific `conda-lock` file before production release and retain it with the installation.
3. Install the package in that environment with `python -m pip install --no-deps .`.
4. Build separate shared references for human GRCh38 and mouse GRCm39. Never mix FASTA, GTF, STAR index, chromosome sizes or PAS atlas releases.
5. Give the server user group read/execute access to software and references. Users need write access only to their own output and scheduler-log directories.
6. Expose the environment with the site module system or a small launcher. Do not modify users' shell startup files.

Example Lmod modulefile logic:

```lua
help([[rna_ends2tracks 0.1.0: QuantSeq REV DGE and APA workflow]])
whatis("Version: 0.1.0")
prepend_path("PATH", "/opt/conda/envs/rna_ends2tracks-0.1.0/bin")
setenv("RNA_ENDS2TRACKS_HOME", "/opt/rna_ends2tracks/0.1.0")
```

PolyAseqTrap and its pinned DeepIP model are deliberately installed as a separate pilot environment. Configure `apa_b.command_template` to call the locally validated adapter. APA-B refuses to run unless `pilot_accepted: true` is explicitly recorded.

## User procedure

Copy the project configuration, samplesheet and one reference-manifest template into the project. Use absolute paths for production FASTQs and references. Validate before submitting compute jobs:

```bash
rna-ends2tracks --config project.yaml --samplesheet samples.csv validate
rna-ends2tracks --config project.yaml --samplesheet samples.csv --dry-run all
rna-ends2tracks --config project.yaml --samplesheet samples.csv all
```

For a scheduler, invoke individual modules in dependency jobs (`preprocess`, then `dge`/`apa-a`/`apa-b`, then `compare` and `report`). Every user runs the same executable; project configuration supplies species and paths.

Use a cooperative project umask such as `002` when a group must edit outputs. The workflow itself never recursively changes permissions and never removes an output tree.
