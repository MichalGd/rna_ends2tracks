# Combined mouse mESC degron QuantSeq run 8 draft

This is one 18-sample mouse project with six conditions: `CTCF_control`,
`CTCF_IAA`, `RAD21_control`, `RAD21_IAA`, `WAPL_control`, and `WAPL_IAA`.
Every condition has three biological replicates. The workflow will generate all
15 pairwise comparisons, including the three matched degron comparisons and all
cross-target comparisons.

## Normalizations and assumptions

- The samples are confirmed as mouse mESC data; the GRCm39/gencode vM31 shared
  reference manifest is selected.
- Libraries are confirmed no-UMI Lexogen QuantSeq REV V2 single-end.
- Every row declares `GRCm39`, matching the selected reference manifest. Each
  biological sample has one library preparation (`technical_replicate_id: T01`)
  and one sequencing lane (`L002`); additional technical preparations or lanes
  would repeat the same `sample_id` in additional rows.
- `description` provides a human-readable label, while filesystem-safe IDs remain
  stable machine identifiers.
- The model is `~ condition`. All samples are from `run8`, so a one-level batch
  term is retained only as metadata and excluded from the statistical design.
- Control and IAA samples are paired biological units within each target. Matching
  subject IDs are recorded as `CTCF_R1` through `CTCF_R3`, `RAD21_R1` through
  `RAD21_R3`, and `WAPL_R1` through `WAPL_R3`.
- Read length is 101 nt, measured directly in representative run-8 FASTQs; the
  table's 100 bp value appears nominal.
- RAD21 replicate IDs follow the replicate column (1, 2, 3), not the inconsistent
  R2/R3/R4 text in `Sample_Name`.
- `WAPL_ctrn` is normalized to `WAPL_control`; lane is recorded as `L002`.
- APA-B remains disabled. Without a PAS atlas, APA-A retains de novo and gene-body
  PCPA analysis but has reduced known-PAS annotation.

Validate on the server with:

```bash
rna-ends2tracks \
  --config mesc_degron.project.yaml \
  --samplesheet mesc_degron.samplesheet.csv \
  validate
```

Validation should report 18 samples, 18 lanes and 15 contrasts. It should mark the
three within-target contrasts as paired and cross-target contrasts as unpaired.

The alpha.5 configuration uses contrast-specific pairing. The three complete
within-target pairs resolve to `~ subject + condition`; the twelve cross-target
contrasts have disjoint subject sets and resolve to `~ condition`. A single global
`~ subject + condition` formula must not be substituted because subjects are nested
within target and that global model is rank-deficient.
