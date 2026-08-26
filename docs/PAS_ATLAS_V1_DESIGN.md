# PAS atlas v1 design for GRCh38 and GRCm39

## Decision

The first production PAS rescue atlases are conservative, assembly-locked administrator assets. They do not define the project PAS universe. Project-supported novel sites remain discoverable condition-blind from QuantSeq data and are never written back into the shared atlas.

| Profile | Target annotation | Core evidence | Rescue evidence |
|---|---|---|---|
| Human | GRCh38 / GENCODE v42 | GENCODE `polyA_site` and PolyA_DB v4.1 Main | PolyA_DB v4.1 Max |
| Mouse | GRCm39 / GENCODE vM31 | GENCODE `polyA_site` and lifted PolyA_DB v4.1 Main | lifted PolyA_DB v4.1 Max |

PolyASite v3 is deferred from v1. It is a discovery-oriented catalog and is not needed by the default `PAS_MASK_RESCUE_TIER=core`. A future extended atlas may add it after an explicit stringency/internal-priming policy is validated; that change requires a new atlas ID and new project output directory.

## Differences from the initial review

The initial review proposed GENCODE, PolyA_DB and PolyASite as a three-source hierarchy. Atlas v1 preserves that architecture but makes these corrections:

1. Mouse is locked to GENCODE vM31 because the installed GTF and STAR index use vM31. M39 is also GRCm39 but must not be mixed into the vM31 profile.
2. Only GENCODE `polyA_site` features are cleavage coordinates. `polyA_signal` is an upstream motif and `pseudo_polyA` is not admitted to the strict core.
3. PolyA_DB v4.1 mouse is distributed on mm10/GRCm38, so Main and Max both require unique, strand-preserving lift-over to GRCm39. The earlier review identified this requirement only for PolyASite mouse.
4. PolyA_DB Max supplies the initial rescue tier. PolyASite is deferred rather than silently applying an unvalidated permissive filter.
5. `novel_supported` is a project-level state, not a shared-atlas class.

## Coordinate and evidence rules

- Raw GTF is 1-based inclusive; normalized source BED is 0-based half-open.
- A GENCODE plus-strand site uses the feature's transcript-direction 3' coordinate; a minus-strand site uses the opposite interval edge.
- PolyA_DB `PAS_ID` positions are normalized from one-based identifiers to one-base BED intervals.
- Duplicate PolyA_DB rows caused by multiple source gene assignments collapse to one genomic PAS. Source identifiers remain in provenance, while final genes are reassigned against the exact target GTF.
- Mouse PolyA_DB sites survive only when lift-over yields exactly one mapping and preserves strand.
- Unknown contigs, out-of-bounds sites and records explicitly marked as internal priming are excluded by the atlas builder.
- Clusters within the configured merge distance prefer PolyA_DB Main coordinates, then GENCODE coordinates, then rescue coordinates.

## Confidence tiers

| Confidence | Rule | Installed tier |
|---|---|---|
| A | GENCODE and PolyA_DB Main occur in the same PAS cluster | core |
| B | GENCODE-only or PolyA_DB Main-only cluster | core |
| C | PolyA_DB Max-only cluster | rescue |

The normal workflow uses `core`. `core_plus_rescue` is an explicit sensitivity choice and changes the PAS-universe signature.

## Reproducible preparation

`scripts/bash/download_pas_atlas_v1_sources.sh` downloads the five fixed official v1 inputs and freezes the snapshot directory. Because downloading or using a UCSC chain indicates acceptance of UCSC terms, it requires the explicit `--accept-ucsc-chain-eula` flag.

`scripts/prepare_pas_sources.py` converts official source snapshots into the strict BED contract consumed by `scripts/build_pas_atlas.py`. It:

- filters GENCODE to `polyA_site`;
- extracts and deduplicates PolyA_DB Main and Max;
- requires a chain and UCSC `liftOver` for foreign mouse coordinates;
- writes deterministic gzip files;
- records raw and normalized SHA-256 checksums, source URLs, releases, coordinate policy and lift-over counts;
- refuses to overwrite an existing prepared-source directory.

The official source profiles are:

- GENCODE v42 human polyA GTF and PolyA_DB v4.1 `HumanPas.zip`;
- GENCODE vM31 mouse polyA GTF and PolyA_DB v4.1 `MousePas.zip`;
- UCSC `mm10ToMm39.over.chain.gz` for mouse coordinate conversion.

The chain file and the `liftOver` executable are separate assets. The preparer accepts an administrator-provided executable through `--lift-over`; it does not modify the immutable workflow environment to install one.

Downloading or using UCSC chain files is subject to the UCSC Genome Browser EULA. An administrator must review and accept those terms; the workflow does not imply acceptance on a user's behalf.

## Release outputs

Each atlas directory contains:

- `core.bed.gz`;
- `rescue.bed.gz`;
- `master.tsv.gz`;
- `provenance.json`;
- `build_report.tsv`;
- `SHA256SUMS`.

Prepared-source directories additionally contain normalized GENCODE/Main/Max BED files, `source_manifest.json`, `preparation_provenance.json` and their checksums. Both source snapshots and final atlas directories are immutable shared assets.

## Promotion gate

Before a workflow release becomes the stable launcher:

1. build human and mouse atlases from immutable official snapshots;
2. verify all checksums and lift-over reports;
3. run one GRCh38 and one GRCm39 metadata/dry-run canary using the installed atlas directories;
4. retain the prior stable launcher as rollback.
