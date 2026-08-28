# APA-B speed-up rationale

## Does reusing C1+C1S compromise APA-B independence?

No. APA-B still detects polyadenylation sites independently. C1+C1S contains only strand-specific, single-nucleotide read endpoints; it does not contain APA-A-called polyadenylation sites.

APA-B independently performs:

- PolyAseqTrap weighted clustering into candidate PAS/PACs;
- its own deterministic summit selection;
- DeepIP internal-priming filtering;
- annotation and support filtering;
- DRIMSeq/stageR differential APA analysis.

APA-B does not reuse APA-A's C2 internal-priming filtering, PAS-atlas rescue decisions, active-PAS catalog, neighbouring-site merging, gene assignments, or statistical results. Those are the operations that define the APA-A site universe.

PolyAseqTrap similarly separates input preparation from PAS clustering, internal-priming assessment, annotation, and quantification. See the [official PolyAseqTrap repository](https://github.com/APAexplorer/PolyAseqTrap).

## Shared preprocessing boundary

The architecture is:

```text
Shared STAR alignment
        |
        v
Raw transcript-oriented endpoints (C1 + C1S)
        |
        +-- APA-A: Mcell2019 internal-priming filtering and PAS calling
        |
        +-- APA-B: PolyAseqTrap clustering and DeepIP filtering
```

In `rna_ends2tracks` terminology:

- C1 contains eligible, unclipped transcript-oriented endpoints;
- C1S contains eligible endpoints whose defining read end is soft-clipped;
- C1+C1S reconstructs the unfiltered C0 endpoint universe used by APA-B;
- C2 is APA-A-specific because it applies internal-priming filtering and rescue logic.

Only C1+C1S may be reused by APA-B. C2 and every later APA-A product are prohibited as APA-B inputs.

This means that APA-A and APA-B are independent PAS callers operating on common validated preprocessing, analogous to two statistical methods operating on the same alignment evidence.

## Why reuse is faster

Without reuse, APA-B reads every complete BAM again and repeats transcript-oriented endpoint extraction. For large QuantSeq projects this can require several hours and substantial shared-storage traffic.

With `APA_B_ENDPOINT_SOURCE="auto"`, APA-B:

1. verifies the successful per-sample exact-end receipt;
2. validates the recorded source files and checksums;
3. verifies the `C0 = C1 + C1S` count invariant;
4. combines C1 and C1S without filtering or deduplication;
5. starts independent PolyAseqTrap/DeepIP processing from the reconstructed raw endpoints.

The optimized and BAM-derived endpoint representations must pass the synthetic coordinate, strand, and count-equivalence pilot before the optimized adapter can be accepted.

## Remaining limitation

Shared endpoint extraction can propagate the same coordinate-extraction software defect into both branches. It does not couple their PAS-calling algorithms, but it reduces engineering independence at the preprocessing boundary.

For a strict preprocessing cross-check, configure:

```bash
APA_B_ENDPOINT_SOURCE="bam"
```

This makes APA-B rescan the BAM files. It should produce the same endpoint evidence and PAS results, but it is slower. A BAM rescan does not make the downstream PolyAseqTrap site-calling algorithm more biologically independent; it provides an additional implementation audit.

## Recommended policy

- Use `APA_B_ENDPOINT_SOURCE="auto"` for routine production analyses.
- Use `APA_B_ENDPOINT_SOURCE="bam"` for initial method validation and occasional audit runs.
- Keep the C1+C1S-versus-BAM synthetic equivalence test mandatory.
- Periodically compare catalogs generated through both endpoint-source modes on representative human and mouse QuantSeq REV data.
- Treat APA-A and APA-B catalogs and statistical results as separate outputs; compare them only in the dedicated downstream concordance module.

Therefore, C1+C1S reuse removes duplicated input processing without replacing or importing APA-A site detection into APA-B.
