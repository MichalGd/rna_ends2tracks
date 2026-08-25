# PolyAseqTrap adapter and REV pilot contract

PolyAseqTrap does not provide a stable QuantSeq REV command interface that this repository can safely guess. A site administrator must pin its source commit, DeepIP model files and dependencies, then provide an adapter command through `apa_b.command_template`.

Available placeholders are `{bam_manifest}`, `{fasta}`, `{gtf}`, `{outdir}`, `{species}` and `{assembly}`. The adapter must preserve exact REV R1 end coordinates, must not extract UMIs, and must not collapse coordinate-identical reads.

It must create these independent files under `{outdir}`:

- `pas_catalog.tsv`: columns `pas_id`, `gene_id`, `chrom`, `start`, `end`, `strand`, `feature_class`, plus PolyAseqTrap/DeepIP evidence columns.
- `pas_counts.tsv`: `pas_id` followed by exactly the biological sample IDs in validated order; values are non-negative raw integers.
- `deepip_audit.tsv`: all candidate IDs, probabilities/classes, model identifier and correction decision, including rejected sites.

No APA-A catalog, identifiers or internal-priming classifications may enter the adapter. After the pilot demonstrates strand/base correctness, count conservation, retention of coordinate-identical reads, genuine/artifact A-rich separation and reproducibility, set `pilot_accepted: true` in the project configuration.
