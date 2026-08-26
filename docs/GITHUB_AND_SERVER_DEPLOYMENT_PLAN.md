# GitHub publication and shared-server deployment plan

> Historical alpha.1 planning record. The current alpha.6 procedure is [server_installation.md](server_installation.md).

## 1. Deployment decision

Publish the current implementation as a pre-release, initially `v0.1.0-alpha.1`, rather than a production release. The Python contracts have passed compilation and smoke testing, but the R modules, Linux bioinformatics tools, synthetic truth set, human/mouse reference integration and PolyAseqTrap REV pilot still require server validation.

Use immutable, side-by-side installations:

```text
/shared/apps/rna_ends2tracks/
├── releases/
│   ├── 0.1.0-alpha.1/       # source from one Git tag
│   └── future-version/
├── envs/
│   ├── 0.1.0-alpha.1/       # independent Conda prefix
│   └── future-version/
├── wheels/
├── modulefiles/
└── logs/

/shared/references/rna_ends2tracks/
├── human/GRCh38/<provider-release>/
└── mouse/GRCm39/<provider-release>/
```

Replace `/shared/apps` with the server's established application root. Software must be read-only to ordinary users. References are independently versioned and are never stored in GitHub.

## 2. Can installation run while analyses are active?

Yes, provided the new installation is isolated from every path used by active jobs.

### Safe during active analyses

- Create or clone a new Git repository in a new directory.
- Download a tagged source archive into a staging directory.
- Build a new Conda environment under a previously unused prefix.
- Build a wheel and run Python/unit tests in the new environment.
- Create a new versioned modulefile without changing the current default.
- Run lightweight validation on a build node or through the scheduler.

### Safe only after considering resource contention

- Conda solving and package extraction can generate substantial shared-filesystem metadata and I/O.
- Container construction can use CPU, storage and network bandwidth.
- Synthetic STAR alignment and R/Bioconductor tests should run as scheduled jobs.
- STAR genome-index generation is CPU-, RAM- and I/O-intensive and should not run interactively beside production jobs.

These operations do not corrupt active analyses when paths are isolated, but they can slow them. Submit high-resource installation checks through the scheduler with explicit limits or wait for a low-load period.

### Unsafe while the affected installation is in use

- Running `conda update`, `pip install`, or `pip uninstall` inside an environment used by active jobs.
- Editing a shared source checkout from which active Python processes or later subprocesses load code.
- Rebuilding or replacing a FASTA, GTF, STAR index, PAS atlas or chromosome-size file in place.
- Changing a shared `current` symlink or unversioned launcher before active workflows have finished.
- Removing old environments, package caches or references.
- Running `conda clean --all` during environment construction or active analyses.

An already running process may have loaded some libraries, but later pipeline steps can start new executables or import additional modules. In-place changes can therefore affect a job hours after it started.

## 3. GitHub repository preparation

### 3.1 Repository boundary

Initialize Git inside the existing `rna_ends2tracks/` directory, not in its parent. This avoids accidentally publishing the reviewed repositories, local analysis documents, datasets or unrelated files.

Proposed metadata:

- Repository: `rna_ends2tracks`
- Default branch: `main`
- Initial visibility: private until the release audit is complete
- Initial tag: `v0.1.0-alpha.1`
- Description: “Portable no-UMI Lexogen QuantSeq REV workflow for DGE, APA and candidate intragenic PCPA analysis in human and mouse.”

### 3.2 Files to add before publication

- Choose and add an explicit `LICENSE`.
- Add `CITATION.cff` with authors, repository URL, version and preferred citation.
- Add `CHANGELOG.md` documenting the alpha scope and known gates.
- Add `CONTRIBUTING.md`, a minimal code of conduct and a security/contact policy if the repository will accept external contributions.
- Add acknowledgements for `rnaseq2tracksP`, `3end-RNAseq-0.1`, Lexogen guidance and the statistical/software methods. Distinguish conceptual influence from copied code.
- Add an output/interface version to every machine-readable contract.
- Decide whether GitHub Issues and Discussions should be enabled.

Do not add FASTQs, BAMs, genome files, annotations, PAS atlases, STAR indexes, model weights, Conda package caches, result directories or credentials. Large biological assets belong in versioned server reference storage or an external archival repository.

### 3.3 Pre-commit release audit

Run from `rna_ends2tracks/`:

```bash
git init -b main
git status --short
git add --all
git diff --cached --check
git diff --cached --stat
```

Before committing:

1. Inspect every staged path.
2. Scan for passwords, tokens, private URLs, email addresses that should not be public and personal filesystem paths.
3. Confirm executable bits on `bin/*.sh`, `scripts/bash/*.sh` and `tests/*.sh`.
4. Normalize text files to LF through `.gitattributes`.
5. Verify that `.gitignore` excludes generated results, caches, environments and build artifacts.
6. Confirm all source links and release metadata.

Then create the initial commit. Creating the remote and pushing should be a separate, explicit action after reviewing the staged tree.

Example remote creation after an owner and visibility have been chosen:

```bash
gh repo create OWNER/rna_ends2tracks --private --source=. --remote=origin --push
```

Do not execute this example until the destination organization/account and visibility are confirmed.

### 3.4 GitHub protections and automation

Configure a repository ruleset for `main`:

- changes enter through pull requests;
- required CI checks must pass;
- force pushes and branch deletion are blocked;
- at least one review is required once there is more than one maintainer;
- tagged releases are created from reviewed commits.

Initial Linux CI should perform:

1. Python 3.10 and 3.11 installation, compilation and unit tests.
2. JSON schema and example YAML/CSV validation.
3. `bash -n` and ShellCheck.
4. R parsing and package-availability checks for DESeq2, DEXSeq, DRIMSeq and stageR.
5. Wheel and source-distribution builds followed by installation into a clean environment.
6. A check that installed package data includes the three R scripts and adapter FASTA.
7. Secret scanning and dependency review.

Add the small synthetic genome/alignment test as the next required CI gate. Full public datasets are too large for ordinary pull-request CI and should run as scheduled release-validation jobs.

### 3.5 Release artifacts

For `v0.1.0-alpha.1`, publish:

- Git source tag and generated source archive;
- Python wheel and source distribution;
- SHA-256 checksum manifest;
- Linux `conda-lock` file for the target architecture;
- example human and mouse manifests without reference assets;
- release notes listing completed tests and unresolved gates.

Do not mark APA-B as validated in the release notes until its pinned PolyAseqTrap/DeepIP pilot has passed.

## 4. Server installation plan

### 4.1 Inventory and decisions

Record before installation:

- Linux distribution and CPU architecture;
- scheduler and permitted build nodes/partitions;
- application and reference roots;
- Conda/Mamba/Micromamba or Apptainer policy;
- existing Lmod/Environment Modules conventions;
- bioinformatics administrator group;
- available RAM/scratch space for STAR indexes;
- whether compute nodes have internet access;
- currently running workflows and the exact environments/references they use.

Recommended first deployment: one immutable Conda prefix created from a Linux lock file, exposed through a versioned modulefile. An Apptainer image can be added after the Conda-based synthetic baseline succeeds.

### 4.2 Build in staging

1. Fetch the exact Git tag into a new staging directory.
2. Verify the Git commit/tag and release checksum manifest.
3. Create a new environment prefix; never activate and modify an existing production prefix.
4. Install the release wheel with `python -m pip install --no-deps`, not editable mode.
5. Verify `rna-ends2tracks --version` and installed R/resource assets.
6. Record the environment package list and executable paths.

Conceptual commands, adjusted for the local Mamba installation:

```bash
micromamba create --yes --prefix /shared/apps/rna_ends2tracks/envs/0.1.0-alpha.1 \
  --file conda-lock-linux-64.yml

/shared/apps/rna_ends2tracks/envs/0.1.0-alpha.1/bin/python -m pip install \
  --no-deps /shared/apps/rna_ends2tracks/wheels/rna_ends2tracks-0.1.0-py3-none-any.whl
```

The exact lock filename and package manager syntax must match the generated lock format.

### 4.3 Reference deployment

Build references in new release-specific directories. For each species:

- record provider, annotation release, assembly and download URLs;
- verify downloaded checksums;
- ensure the FASTA, GTF, PAS atlas and chromosome naming agree;
- build `.fai`, chromosome sizes and STAR index;
- record STAR version and index-generation command;
- make the completed directory read-only;
- validate it with `rna-ends2tracks validate` using a small species-specific project.

Before building a new STAR index, audit any existing index. Reuse is preferred when its required files are complete, contig names and lengths exactly match the chosen FASTA, FASTA/GTF provenance is documented, STAR compatibility and `sjdbOverhang` are acceptable, and the directory is immutable to users. Concurrent read-only index use is correctness-safe, although simultaneous index loading can add storage I/O. A mismatched, composite, incomplete, mutable or provenance-incompatible index must not be reused.

Do not share one unversioned `GRCh38` or `GRCm39` directory among releases if files can later be replaced. Publish a new directory and update only new project manifests.

### 4.4 Versioned modulefile

Create `rna_ends2tracks/0.1.0-alpha.1`. It should prepend only the new environment's `bin` directory and set a read-only workflow-home variable. It must not alter Python library paths globally or select a reference automatically.

Users should explicitly load the alpha version:

```bash
module load rna_ends2tracks/0.1.0-alpha.1
rna-ends2tracks --version
```

Do not make the alpha module the site default. After release validation, promote a stable version through the module system without deleting the alpha or previous stable version.

### 4.5 Permissions

- Application source, wheels, environments and references: owned by the administrator/service account and readable/executable by the bioinformatics user group.
- Modulefiles: readable by all intended users.
- Project results: written to user/group project storage, never under the application installation.
- Reference and application trees: ordinary users must not have write permission.
- Shared projects may use an agreed group and umask, but the installer should not recursively change unrelated project permissions.

### 4.6 Canary validation

Run as an ordinary server user, not as the installer:

1. CLI version/help and metadata-only validation.
2. Python unit tests and Bash syntax checks.
3. R package loading and parsing of all three R modules.
4. Tiny synthetic human-like and mouse-like reference tests.
5. Lexogen REV orientation and exact plus/minus coordinate truth tests.
6. Duplicate-flagged coordinate-identical read retention.
7. Two-lane-to-one-sample merge and all-pairs contrast generation.
8. featureCounts/DESeq2 and APA-A DEXSeq known-effect tests.
9. BigWig generation and coordinate inspection.
10. Receipt-based restart and `--force-module` behavior.

Run one small real QuantSeq REV project after the synthetic suite. Compare mapping, gene counts, known PAS overlap, motif enrichment and representative loci with the prior workflow.

APA-B has an additional, independent canary: pin PolyAseqTrap and DeepIP assets, satisfy the adapter contract, prove no hidden UMI/deduplication behavior and meet the REV pilot criteria before changing `pilot_accepted`.

## 5. Rollout and rollback

### Rollout

1. Keep the release private and install it as an alpha module.
2. Grant access to one or two canary users.
3. Run human and mouse synthetic tests, then one small real project per available species.
4. Record defects and fixes through pull requests; create a new tag for every installed change rather than editing the installed release.
5. After acceptance, tag a release candidate, rebuild from scratch and repeat the canary.
6. Promote a stable module only after all primary gates pass.

### Rollback

Rollback consists of loading the previous module version and using its unchanged reference manifest. No environment needs to be repaired in place.

Keep at least the previous validated software environment, wheel, source tag, lock file and reference manifest until all jobs launched with it have completed and their retention period has expired.

## 6. Acceptance checklist

- [ ] GitHub owner, visibility and license approved.
- [ ] Staged tree contains only the intended `rna_ends2tracks` repository.
- [ ] No secrets, datasets, references, results or personal paths are committed.
- [ ] Linux CI passes Python, Bash, R, package-data and build tests.
- [ ] Tagged wheel installs into a clean locked environment.
- [ ] Versioned module works for an unprivileged user.
- [ ] Application and reference directories are immutable to users.
- [ ] Human GRCh38 synthetic/canary validation passes.
- [ ] Mouse GRCm39 synthetic/canary validation passes.
- [ ] Active analyses continue using their original environment and references.
- [ ] Restart and rollback tests pass.
- [ ] APA-A biological acceptance tests pass.
- [ ] APA-B remains disabled, or its independent REV pilot is documented as passed.
- [ ] Release notes distinguish candidate PCPA from proven premature transcription termination.

## 7. Authoritative operational references

- GitHub: [Creating a new repository](https://docs.github.com/en/repositories/creating-and-managing-repositories/creating-a-new-repository)
- GitHub: [About repository rulesets](https://docs.github.com/en/repositories/configuring-branches-and-merges-in-your-repository/managing-rulesets/about-rulesets)
- GitHub: [Managing releases](https://docs.github.com/en/repositories/releasing-projects-on-github/managing-releases-in-a-repository)
- Conda: [Managing environments](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html)
- Apptainer: [SIF image format](https://apptainer.org/docs/user/latest/sif.html)
