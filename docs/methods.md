# Methods and count universes

## Preprocessing and alignment

The validated profiles are Lexogen QuantSeq REV V2 single-end and paired-end libraries without UMIs. Raw reads are optionally screened against the site-configured FastQ Screen species/contaminant databases, raw and trimmed mates are inspected with FastQC, and all QC is summarized by MultiQC. BBDuk uses the Lexogen-oriented adapter/poly(A/T) reference with `k=13 ktrim=r useshortkmers=t mink=5 qtrim=r trimq=10 minlength=20` by default. Paired inputs use synchronized `in1/in2` and `out1/out2`. A second synchronized BBDuk pass applies `ftl=12 skipr1=t`, removing Lexogen's random-primer-derived first 12 bases from R2 only; `PE_R2_TRIM_5P` makes this explicit and configurable. All values that are safe to vary are exposed in `config.conf`.

STAR maps reads to the assembly selected in each samplesheet row. Human and mouse samples may coexist, but each lane uses its matching reference and every downstream catalog/statistical model remains genome-specific. The workflow audits STAR `ReadsPerGene.out.tab` and expects the reverse-stranded QuantSeq REV profile. Existing STAR indices are accepted after structural validation; `sjdbOverhang` differences generate a review warning because a 150-overhang index can still map shorter 101-nt reads.

After alignment, RSeQC independently evaluates the final C0 BAMs with `infer_experiment.py`, `read_distribution.py`, and `geneBody_coverage.py`. The reference is an assembly-matched BED12 supplied in `config.conf` or deterministically generated from the exact selected GTF. Per-sample work is bounded by `RSEQC_PARALLEL_JOBS`; its summaries, combined gene-body plot, source files, receipt, and dedicated MultiQC report are retained under `01_qc/rseqc`. QuantSeq REV is expected to show strong transcript 3-prime enrichment, so conventional whole-transcript uniformity is not the acceptance criterion. See [RSeQC QC](rseqc.md).

The alignment BAM contains mapped primary `NH=1` alignments. SAM flags 0x4, 0x100 and 0x800 are excluded. In SE it contains one alignment per mapped read; in PE it retains both mapped mates. Duplicate flags are retained because no UMI establishes molecular identity and real 3-prime molecules often share a cleavage coordinate.

## Exact transcript ends

QuantSeq REV Read 1 is antisense to the source transcript. A reverse genomic R1 alignment represents a transcript-plus molecule and its transcript 3′ coordinate is the rightmost aligned reference base (`reference_end-1`). A forward R1 alignment represents transcript-minus and uses `reference_start`. For paired libraries, both mates contribute their aligned blocks to conventional all-read coverage, which is normalized per mapped pair by counting R1 records. Split CIGAR blocks are honored so introns are not painted as covered sequence. R2 is explicitly excluded from the APA-A/APA-A2 C1/C1S coordinate stream and from APA-B cleavage-coordinate extraction. These tracks show aligned read blocks rather than an inferred insert span.

The defining CIGAR side is checked for soft/hard clipping. Unclipped coordinates form C1; uncertain clipped records form C1S. Mapping-class, clipping and duplicate audits are retained. The end-analysis C0 universe means eligible end-defining molecules (R1 fragments in PE), not the number of mate alignments in a PE BAM. Required identity: `C0 = C1 + C1S`.

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

C5 reverse-stranded featureCounts exon counts are diagnostic only. Single-end projects count reads; paired-end projects use `-p --countReadPairs` and count fragments. The workflow reports per-sample log-count correlation and genes whose C4/C5 CPM ratio differs by at least fourfold. It never silently substitutes C5 for C4.

## APA-A and corrected APA-A2

Legacy APA-A is preserved unchanged. Genes with at least two eligible active PAS are tested by DEXSeq using raw C3 counts; ambiguous sites are excluded. Its historical effect and direction layer derives PAU from DEXSeq-normalized counts and classifies a two-site comparison. If at least two sites are significant, the two lowest adjusted p-values are chosen, then absolute delta-PAU, mean normalized count and coordinate break ties. If exactly one is significant, its comparator is the nonzero other PAS with highest mean normalized usage, then pooled raw count and coordinate. Proximal/distal order follows transcript direction and its zero-safe ratio direction uses a cross product.

APA-A2 independently reruns the same DEXSeq hypothesis test on the shared, condition-blind C3 catalog. Its effect layer calculates within-gene PAS usage directly from raw C3 counts separately in every sample. Zero-total gene/sample observations remain `NA`. Unpaired delta-PAU is the difference between condition means; paired delta-PAU is calculated within each complete subject and averaged with equal pair weight. A primary site requires both site FDR and `MIN_ABS_DELTA_PAU`; a primary gene requires gene-level FDR and at least one primary site. Direction uses the PAU-weighted mean transcript-coordinate shift over the complete selected site set, with transcript-minus coordinates reversed. Every contrast publishes PAU-sum and within-gene effect-conservation audits. See [APA-A2 corrected analysis](12_apa_a2_corrected.md).

## Pairwise designs

All within-genome condition pairs meeting the configured biological-replicate minimum are generated. Complete one-to-one shared subjects use `~ subject + condition`; disjoint subjects use `~ condition`. Partial overlap/incomplete pairs fail by default. Technical libraries and lanes collapse before replication is counted.
