# rna_ends2tracks 2.0

The release name is **2.0**, the Python/CLI version is **2.0.0**, and the canonical Git tag is **v2.0.0**. Tag the commit containing this version after the updated code has been merged and checked. Do not attach the tag to an older commit whose package still identifies itself as an alpha.

This names the alpha13.post1 implementation already installed and verified on the shared server. It retains PAS-independent extended-exon counting, depth-aware expression filtering, separate stringent/relaxed DE and ORA, shared signed GSEA, BAM-input reanalysis and distance-clustered VST heatmaps. The quantification/statistical algorithms and dependency requirements are unchanged in the version promotion. See [the expression method guide](gene_expression_alpha13.md) and [the changelog](../CHANGELOG.md).

## Shared-server layout

| Component | Location |
|---|---|
| Existing dependency environment, shared by alpha12 and 2.0 | `/opt/conda_envs/rna_ends2tracks-0.1.0a12.post1` |
| Small 2.0 workflow-code layer | `/opt/rna_ends2tracks/2.0.0` |
| Versioned workflow launcher | `/opt/conda_envs/bin/rna-ends2tracks-2.0.0` |
| Shared default command | `rna-ends2tracks` |

The application layer uses the existing Python interpreter and packages; R and command-line tools come from the same dependency environment. A local-wheel install needs no dependency download, Conda creation or full-environment copy. Keep the dependency environment at its existing path. Versioned alpha12/alpha13.post1 commands remain available for previous analyses. Each user's project controls its output and temporary storage; run11's `/dysk3` choice is independent of the software installation on the SSD.

After deploying the 2.0 wheel and switching the shared launcher, `rna-ends2tracks --version` must print `2.0.0`. A symlink alias alone does not change package metadata, reported versions or provenance.

## Results and compatibility

2.0 accepts a successful alpha13.post1 receipt only if the full original signature and output-file integrity checks also pass. Earlier alpha receipt families are not declared equivalent to the new expression implementation. Changes to data, references, configuration or scripts can still invalidate individual cached stages. Read-only BAM reuse into a separate expression project remains the intended way to redo the older run11 analyses. Version promotion does not validate the full-depth biological PMK pilot; that remains a separate analytical check.

## Publication sequence

1. Build/check the `2.0.0` wheel and reviewed source snapshot.
2. Install the small 2.0 application layer over the existing shared dependency environment; confirm the version and ordinary-user access.
3. Publish the reviewed source changes to a branch and merge after checks pass.
4. Attach annotated tag `v2.0.0` to that merged commit and push the tag. Verify it identifies the reviewed source tree.
5. Create a GitHub release named `rna_ends2tracks 2.0` and attach the wheel, source snapshot, checksum list and release notes.

This document defines names and sequencing; it does not assert that the server or remote repository has already been updated to 2.0.
