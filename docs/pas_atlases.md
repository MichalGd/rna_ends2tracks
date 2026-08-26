# Human and mouse PAS atlases

The production v1 evidence policy and its differences from the initial method review are defined in [PAS_ATLAS_V1_DESIGN.md](PAS_ATLAS_V1_DESIGN.md).

PAS atlases are versioned administrator-built assets used only to rescue plausible true cleavage sites from the internal-priming mask. They do not define the project's active PAS universe, and novel project-supported sites remain eligible.

Initial profiles are:

| Species | Assembly | Annotation | Atlas ID |
|---|---|---|---|
| human | GRCh38 | GENCODE v42 | `GRCh38_gencode_v42_pas_atlas_v1` |
| mouse | GRCm39 | GENCODE vM31 | `GRCm39_gencode_vM31_pas_atlas_v1` |

Each atlas directory contains `core.bed.gz`, `rescue.bed.gz`, `master.tsv.gz`, `provenance.json`, `build_report.tsv` and `SHA256SUMS`. Core A is concordant GENCODE/PolyA_DB Main evidence. Core B is assembly-native GENCODE or PolyA_DB Main evidence. Rescue C is PolyA_DB Max evidence absent from the core. Multi-gene status is retained in the master table.

## Download immutable official snapshots

```bash
bash scripts/bash/download_pas_atlas_v1_sources.sh \
  --output /snapshots/PAS_atlas_v1_raw_YYYYMMDD \
  --accept-ucsc-chain-eula
```

The acceptance flag records an administrator decision; read the linked UCSC EULA before using it. The download script refuses to replace an existing snapshot and makes a completed snapshot read-only.

Mouse normalization also requires a UCSC `liftOver` executable on `PATH`, or an explicit executable path supplied with `--lift-over`.

## Normalize sources

Mouse:

```bash
python scripts/prepare_pas_sources.py \
  --species mouse --assembly GRCm39 --annotation-release GENCODE_vM31 \
  --gencode-polya-gtf /snapshots/PAS_atlas_v1_raw_YYYYMMDD/gencode.vM31.polyAs.gtf.gz \
  --polyadb-zip /snapshots/PAS_atlas_v1_raw_YYYYMMDD/MousePas.zip \
  --polyadb-assembly mm10 \
  --liftover-chain /snapshots/PAS_atlas_v1_raw_YYYYMMDD/mm10ToMm39.over.chain.gz \
  --download-date YYYY-MM-DD \
  --output /snapshots/prepared_mouse_GRCm39_v1
```

Human:

```bash
python scripts/prepare_pas_sources.py \
  --species human --assembly GRCh38 --annotation-release GENCODE_v42 \
  --gencode-polya-gtf /snapshots/PAS_atlas_v1_raw_YYYYMMDD/gencode.v42.polyAs.gtf.gz \
  --polyadb-zip /snapshots/PAS_atlas_v1_raw_YYYYMMDD/HumanPas.zip \
  --polyadb-assembly hg38 \
  --download-date YYYY-MM-DD \
  --output /snapshots/prepared_human_GRCh38_v1
```

## Build an atlas

Mouse example:

```bash
python scripts/build_pas_atlas.py \
  --species mouse --assembly GRCm39 --annotation-release GENCODE_vM31 \
  --atlas-id GRCm39_gencode_vM31_pas_atlas_v1 \
  --source-manifest /snapshots/prepared_mouse_GRCm39_v1/source_manifest.json \
  --gtf /refs/gencode.vM31.gtf --chrom-sizes /refs/GRCm39.chrom.sizes \
  --gencode /snapshots/prepared_mouse_GRCm39_v1/gencode_polyA.GRCm39.bed.gz \
  --polyadb-main /snapshots/prepared_mouse_GRCm39_v1/polyadb_v4.1_main.GRCm39.bed.gz \
  --polyadb-max /snapshots/prepared_mouse_GRCm39_v1/polyadb_v4.1_max.GRCm39.bed.gz \
  --output /shared/references/mouse/GRCm39/GRCm39_gencode_vM31_pas_atlas_v1
```

Use the corresponding human prepared paths, GRCh38 GTF/chromosome sizes and atlas ID for the human build.

Prepared inputs are BED6-plus with 0-based half-open coordinates. Mouse PolyA_DB coordinates require UCSC `liftOver`; only one-to-one, strand-preserving mappings survive. Out-of-bounds, unknown-contig and internal-priming-flagged sites are rejected. Every retained cluster is reannotated against the exact target GTF. PolyASite is deferred from v1 and remains an optional input for a separately versioned future atlas.

The generated `source_manifest.json` contains one complete entry for every normalized GENCODE and PolyA_DB input, including its checksum and raw-source provenance. `preparation_provenance.json` separately records the raw archives and the mouse lift-over chain. Every checksum must match before building.

Changing an atlas, rescue tier, annotation or assembly changes the active-PAS signature. Use a new output directory for the resulting PAS universe. Source snapshots and completed atlas directories are immutable shared resources.
