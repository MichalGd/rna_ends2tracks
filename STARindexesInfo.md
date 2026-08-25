# STAR indexes and single-end QuantSeq REV support

## One STAR index can serve several read lengths

Separate STAR indexes are not required for 50 bp, 100 bp and 150 bp reads when the assembly and annotation are unchanged.

- STAR recommends `sjdbOverhang = max(read length) - 1` when datasets have different read lengths.
- The existing human GRCh38 and mouse GRCm39 indexes use `sjdbOverhang=149`, so they are the shared candidates for reads up to 150 bp.
- Adapter and quality trimming naturally creates shorter reads and does not require another index.
- A new index should be considered for reads longer than 150 bp or when the FASTA/GTF changes, not merely because reads are shorter.
- `rna_ends2tracks` records a review warning when a project's declared maximum read length does not exactly match the index overhang. It does not silently rebuild or switch indexes.

STAR's official manual states that the ideal value is read length minus one, that varying-length datasets should use the maximum read length minus one, and that the default often performs as well as the ideal value: <https://github.com/alexdobin/STAR/blob/master/extras/doc-latex/STARmanual.tex>.

## Single-end contract

`rna_ends2tracks` is intentionally restricted to no-UMI, single-end Lexogen QuantSeq REV data in the validated alpha profiles.

- `library_layout` must be `SE`.
- `fastq_r2` must be empty.
- Enabled profiles are `quantseq_rev_v1_se` and `quantseq_rev_v2_se`.
- STAR receives one trimmed R1 FASTQ per lane.
- Technical lanes belonging to one sample are merged after alignment.
- Paired-end input is rejected until a separate REV R2 and end-coordinate pilot is completed.

## Audited server indexes

On `biolserv`, the existing GRCh38/GENCODE v42 and GRCm39/GENCODE vM31 STAR indexes passed exact FASTA contig-name and contig-length comparison. Both record `sjdbOverhang=149` and the intended GTF path. The indexed human FASTA copy was also byte-identical to the FASTA path recorded during index construction.

These indexes can therefore be reused for the alpha canary. The release manifests must retain the audited FASTA, GTF and index pairing, record checksums and use chromosome-size files derived from the complete matching FASTA indexes. The older `hs38n.chrom.sizes` and `mm39n.chrom.sizes` subset files are not used by `rna_ends2tracks`.

## Alpha.1 preprocessing defect

Review during installation found that `v0.1.0-alpha.1` attempted to read `ReadsPerGene.out.tab` before STAR ran. The orientation check must run after STAR creates that file. The corrected order and a non-dry-run regression test are included in `v0.1.0-alpha.2`; the immutable alpha.1 tag must not be installed for data analysis.
