# Human and mouse PAS atlases

PAS atlases are versioned administrator-built assets used only to rescue plausible true cleavage sites from the internal-priming mask. They do not define the project’s active PAS universe and novel project-supported sites remain eligible.

Initial profiles are:

| Species | Assembly | Annotation | Atlas ID |
|---|---|---|---|
| human | GRCh38 | GENCODE v42 | `GRCh38_gencode_v42_pas_atlas_v1` |
| mouse | GRCm39 | GENCODE vM31 | `GRCm39_gencode_vM31_pas_atlas_v1` |

Each atlas directory contains `core.bed.gz`, `rescue.bed.gz`, `master.tsv.gz`, `provenance.json`, `build_report.tsv` and `SHA256SUMS`. Core A is concordant GENCODE/PolyA_DB Main evidence. Core B is assembly-native GENCODE or PolyA_DB Main evidence. Rescue C is filtered PolyA_DB Max/PolyASite evidence. Multi-gene status is retained in the master table.

Build from immutable local source snapshots:

```bash
python scripts/build_pas_atlas.py \
  --species mouse --assembly GRCm39 --annotation-release GENCODE_vM31 \
  --atlas-id GRCm39_gencode_vM31_pas_atlas_v1 \
  --source-manifest /snapshots/PAS_source_manifest.json \
  --gtf /refs/gencode.vM31.gtf --chrom-sizes /refs/GRCm39.chrom.sizes \
  --gencode /snapshots/gencode_polyA.GRCm39.bed.gz \
  --polyadb-main /snapshots/polyadb_v4_main_mouse.tsv.gz \
  --polyadb-max /snapshots/polyadb_v4_max_mouse.tsv.gz \
  --polyasite /snapshots/polyasite_v3_mouse.GRCm38.bed.gz \
  --polyasite-assembly GRCm38 --liftover-chain /refs/mm10ToMm39.over.chain.gz \
  --output /shared/references/mouse/GRCm39/GRCm39_gencode_vM31_pas_atlas_v1
```

Inputs must be BED6 or headered TSV with 0-based half-open `chrom/start/end/strand`; optional fields are `pas_id`, `gene_id` and `internal_priming_flag`. Foreign mouse PolyASite coordinates require UCSC `liftOver`; only one-to-one, strand-preserving mappings survive. Out-of-bounds, unknown-contig and internal-priming-flagged sites are rejected. Every retained cluster is reannotated against the exact target GTF.

The source manifest is JSON: `{"sources":[{"file":"snapshot.bed.gz","url":"https://…","release":"…","download_date":"YYYY-MM-DD","license":"…","sha256":"…"}]}`. It must contain one complete entry for every GENCODE/PolyA_DB/PolyASite/chain snapshot and every checksum must match before building.

Changing an atlas, rescue tier, annotation or assembly changes the active-PAS signature. Use a new output directory for the resulting PAS universe. Source URLs, release labels and licenses should be stored beside the immutable snapshots; the builder records paths and SHA-256 checksums but does not download data during a user run.
