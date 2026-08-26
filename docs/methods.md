# Methods and count universes

## Preprocessing and alignment

The validated profile is Lexogen QuantSeq REV V2 single-end Read 1 without UMIs. Raw and trimmed reads are inspected with FastQC and summarized by MultiQC. BBDuk uses the Lexogen-oriented adapter/poly(A/T) reference with `k=13 ktrim=r useshortkmers=t mink=5 qtrim=r trimq=10 minlength=20` by default. All values that are safe to vary are exposed in `config.conf`.

STAR maps reads to the assembly selected in each samplesheet row. Human and mouse samples may coexist, but each lane uses its matching reference and every downstream catalog/statistical model remains genome-specific. The workflow audits STAR `ReadsPerGene.out.tab` and expects the reverse-stranded QuantSeq REV profile. Existing STAR indices are accepted after structural validation; `sjdbOverhang` differences generate a review warning because a 150-overhang index can still map shorter 101-nt reads.

The C0 statistical BAM contains mapped primary `NH=1` alignments. SAM flags 0x4, 0x100 and 0x800 are excluded. Duplicate flags are retained because no UMI establishes molecular identity and real 3′ molecules often share a cleavage coordinate.

## Exact transcript ends

QuantSeq REV reads are antisense to the source transcript. A reverse genomic alignment represents a transcript-plus molecule and its transcript 3′ coordinate is the rightmost aligned reference base (`reference_end-1`). A forward alignment represents transcript-minus and uses `reference_start`.

The defining CIGAR side is checked for soft/hard clipping. Unclipped coordinates form C1; uncertain clipped records form C1S. Mapping-class, clipping and duplicate audits are retained. Required identity: `C0 = C1 + C1S`.

## Internal-priming mask

For transcript-plus, a C1 coordinate is masked when it belongs to a genomic window containing six consecutive A bases or at least seven A bases in ten nucleotides. Transcript-minus uses T. A masked coordinate is rescued when it lies in the conservative 20-nt interval around an assembly-matched GENCODE transcript end or installed core PAS-atlas site. `core_plus_rescue` adds the atlas rescue tier and changes the PAS-universe signature.

Mask-pass/rescued records form C2; rejected records form C2R. Required identity: `C1 = C2 + C2R`.

## Mcell2019 PAS discovery

For each sample, C2 is scaled once across both transcript strands:

```text
C2_CPM(position) = 1,000,000 × C2(position) / total sample C2
```

Scaled samples are summed separately by strand within each genome, without using conditions. For each strand, the workflow sums signal in every 30-nt half-open window, retains windows with a sum strictly greater than 30, merges overlapping or adjacent retained windows, and chooses the highest single-base pooled signal in each component. Equal maxima choose the lowest genomic coordinate, reproducing the sorted legacy `[1]` rule on both strands.

For a one-base summit `[p,p+1)`, recentering produces `[p-14,p+16)`. This is the exact 0-based representation of the even-width center convention. It is clipped at chromosome boundaries. The recentered intervals are merged and recentered once more. Final intervals must not overlap.

Raw C2 coordinates are counted once in those intervals to form C3. The catalog stores summit and interval coordinates separately. Required checks are `sum(C3) <= sum(C2)` and at most one C3 assignment per C2 end.

## Gene assignment and premature termination candidates

Genes are extended 6 kb in transcript 3′ direction. Same-strand PAS are classified as terminal exon, other exon, intron or downstream extension. Opposite-strand gene-body sites are antisense; remaining sites are intergenic. Multi-gene overlaps are marked ambiguous and excluded from statistical matrices, while remaining visible in catalogs and tracks.

Unique intronic/non-terminal-exonic PAS in genes that also possess a terminal-exon PAS form the PCPA candidate catalog. Significant candidate PCPA requires DEXSeq-adjusted `p <= FDR` and `|delta_PAU| >= MIN_ABS_DELTA_PAU`. These are candidates consistent with premature cleavage/polyadenylation and premature transcription termination, not direct proof of polymerase termination.

## Gene expression

C4 sums C3 over uniquely assigned PAS per gene and is the primary raw-integer DESeq2 matrix. Genome-global C4 `poscounts` size factors produce visualization scales. Each pairwise contrast is fitted independently using its resolved paired or unpaired design.

C5 reverse-stranded featureCounts exon counts are diagnostic only. The workflow reports per-sample log-count correlation and genes whose C4/C5 CPM ratio differs by at least fourfold. It never silently substitutes C5 for C4.

## APA and shift direction

Genes with at least two eligible active PAS are tested by DEXSeq using raw C3 counts. Ambiguous sites are excluded. If at least two sites are significant, the two lowest adjusted p-values are chosen, then absolute delta-PAU, mean normalized count and coordinate break ties. If exactly one is significant, its comparator is the nonzero other PAS with highest mean normalized usage, then pooled raw count and coordinate.

Proximal/distal order follows transcript direction. Zero-safe direction uses the cross product instead of adding a pseudocount. Both zero is `NA`, nonzero divided by zero is `Inf`, and a library with no classifiable proximal/distal use is reported as not classifiable.

## Pairwise designs

All within-genome condition pairs meeting the configured biological-replicate minimum are generated. Complete one-to-one shared subjects use `~ subject + condition`; disjoint subjects use `~ condition`. Partial overlap/incomplete pairs fail by default. Technical libraries and lanes collapse before replication is counted.
