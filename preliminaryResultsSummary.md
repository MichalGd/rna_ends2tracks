# Preliminary results summary: `mesc_degron_run8`

## Scope and status

This document summarizes the first completed `rna_ends2tracks` analysis of 18 mouse Lexogen QuantSeq REV V2 libraries:

- workflow release: `v0.1.0-alpha.9.post2`;
- genome and annotation: GRCm39 / GENCODE vM31;
- design: CTCF, RAD21 and WAPL control/IAA pairs, three biological pairs per target;
- UMI processing: disabled;
- primary expression universe: C4 active-PAS gene sums;
- primary APA method: APA-A Mcell2019-style active-PAS analysis;
- APA-B: disabled for this run.

The workflow completed successfully. All required stage receipts are present, all 15 configured DGE and APA-A comparisons finished, 504 BigWigs were generated, the basic alpha.9 HTML report was written, and success-gated cleanup removed approximately 111 GiB of dispensable intermediates. No error was found in the final `alpha.9.post2` run log.

These conclusions are preliminary. They establish technical performance and identify patterns worth investigating, but they do not replace inspection of MultiQC, PCA, effect sizes, replicate-level behavior, genomic tracks, or orthogonal validation of degron efficiency.

## Plain-language explanation of C2

`rna_ends2tracks` uses the labels C0–C5 for successive count universes. They are workflow-specific definitions, not universal RNA-seq terminology. The letter C means **count universe**: each successive label identifies a precisely defined set of reads, end coordinates or aggregated counts.

| Universe | Meaning | Main use |
|---|---|---|
| C0 | Mapped, primary, `NH=1` alignments; duplicate flags are retained | Final BAMs and all-read coverage tracks |
| C1 | Exact transcript-end coordinates from eligible C0 alignments whose end-defining coordinate is not obscured by clipping | Raw exact-end signal |
| C1S | C0 alignments clipped at the side needed to define the exact transcript end | Separate uncertainty/QC output; excluded from PAS calling |
| **C2** | **C1 exact ends that pass the internal-priming filter, plus filtered sites rescued by assembly-matched transcript-end/PAS evidence** | **Condition-blind active-PAS discovery and filtered-end tracks** |
| C2R | C1 ends classified as probable internal-priming artifacts and not rescued | Diagnostic reject counts/tracks; excluded from PAS discovery and statistics |
| C3 | C2 ends assigned at most once to the final non-overlapping project-active PAS intervals | PAS usage matrix and APA analysis |
| C4 | C3 counts summed per gene over uniquely assigned PAS | Primary gene-expression matrix and normalization factors |
| C5 | Conventional reverse-stranded exon-overlap featureCounts from C0 BAMs | Diagnostic comparison with C4 only |

### What C2 represents biologically

QuantSeq REV reads are intended to report transcript 3′ ends, but oligo-dT priming can also occur at internal A-rich genomic sequence. A mapped read end is therefore not automatically accepted as a cleavage/polyadenylation site.

For each sample, the workflow first obtains an exact one-base end coordinate (C1). It then applies a strand-aware internal-priming mask around that coordinate. An end becomes C2 when it either:

1. does not satisfy the configured A-rich internal-priming rule; or
2. is A-rich but lies close enough to trusted, assembly-matched transcript-end or PAS-atlas evidence to be rescued.

An A-rich end without rescue evidence becomes C2R instead. Consequently, C2 is the workflow's filtered set of plausible 3′-end observations. It is narrower than all mapped reads and is the direct input to project-wide, condition-blind active-PAS discovery. C2 is not yet a PAS catalog: C2 coordinates must still be pooled, clustered and assigned to final non-overlapping active-PAS intervals to become C3.

## Sequencing depth and mapping

The counts below are reads entering STAR after trimming, not the original raw FASTQ read counts. `Unique C0 reads` are the mapped unique-primary counts retained in the final C0 BAMs.

| Sample | STAR input reads | Unique C0 reads | Unique mapping |
|---|---:|---:|---:|
| CTCF control R1 | 37.69 M | 31.13 M | 82.61% |
| CTCF control R2 | 37.11 M | 30.83 M | 83.07% |
| CTCF control R3 | 32.13 M | 27.12 M | 84.41% |
| CTCF IAA R1 | 29.03 M | 22.99 M | 79.21% |
| CTCF IAA R2 | 21.61 M | 16.93 M | 78.34% |
| CTCF IAA R3 | 20.13 M | 16.71 M | 82.98% |
| RAD21 control R1 | 30.77 M | 23.81 M | 77.40% |
| RAD21 control R2 | 33.23 M | 25.29 M | 76.10% |
| RAD21 control R3 | 34.53 M | 26.55 M | 76.89% |
| RAD21 IAA R1 | 36.71 M | 26.92 M | 73.35% |
| RAD21 IAA R2 | 22.25 M | 16.47 M | 74.03% |
| RAD21 IAA R3 | 28.96 M | 22.43 M | 77.45% |
| WAPL control R1 | 34.51 M | 27.22 M | 78.88% |
| WAPL control R2 | 30.31 M | 23.67 M | 78.09% |
| WAPL control R3 | 33.83 M | 26.91 M | 79.53% |
| WAPL IAA R1 | 26.98 M | 21.35 M | 79.13% |
| WAPL IAA R2 | 27.66 M | 20.87 M | 75.47% |
| WAPL IAA R3 | 14.37 M | 11.37 M | 79.17% |

Aggregate observations:

- 531.8 million post-trimming STAR-input reads;
- mean 29.5 million and median 30.5 million reads per sample;
- range 14.37–37.69 million reads per sample;
- 418.6 million unique C0 reads in total;
- mean 23.3 million unique C0 reads per sample;
- mean unique mapping 78.7%, range 73.35–84.41%.

This is substantial depth for a 3′ RNA-seq experiment, with no outright mapping failure.

## Strand specificity

All 18 samples passed the QuantSeq REV orientation check. The reverse-compatible fraction was:

- mean: 93.07%;
- range: 92.47–93.40%.

The high and unusually consistent reverse-compatible fraction is strong evidence that the libraries behave as strand-specific Lexogen QuantSeq REV libraries. There is no obvious orientation-swapped or effectively unstranded sample.

Strand specificity and good mapping do not by themselves prove correct 3′-end enrichment. That conclusion also depends on C1-to-C2 retention, known/terminal PAS localization, internal-priming rejection, C3 assignment and browser inspection.

## Filtered ends and active PAS

The completed active-PAS audit reported:

| Metric | Result |
|---|---:|
| C2 filtered exact ends | 242,208,168 |
| C3 ends assigned to active PAS | 221,966,192 |
| C3 as a percentage of C2 | 91.64% |
| Project-active PAS | 31,407 |
| Ambiguous PAS | 3,368 (10.72%) |
| Candidate PCPA sites in the catalog universe | 5,944 (18.93% of active PAS) |

Approximately 58% of the aggregate unique-primary C0 alignments remained in C2. This combines exclusion of end-defining clipped reads and rejection of probable internal priming; the per-sample C1/C1S/C2/C2R funnel is still needed to separate those causes.

The assignment of 91.64% of C2 signal to C3 is reassuring: most filtered exact-end signal is represented by the final project-active PAS universe. The audit also confirmed that active-PAS intervals are non-overlapping, C3 does not exceed C2, and each C2 end is assigned no more than once.

The 5,944 PCPA sites are the candidate intragenic PCPA **universe**, not 5,944 significant differential premature-termination events. Significant contrast-specific candidates are reported separately in `06_apa_a_mcell2019/GRCm39/candidate_pcpa.tsv`. QuantSeq evidence supports candidate premature cleavage/polyadenylation but does not by itself prove loss of downstream engaged RNA polymerase.

## Preliminary DGE and APA results

The three paired within-target comparisons are the primary treatment contrasts.

| Primary paired comparison | Significant DGE genes | APA-A tested sites | Significant APA-A sites | Significant fraction |
|---|---:|---:|---:|---:|
| CTCF IAA vs CTCF control | 1,981 | 21,620 | 6,795 | 31.43% |
| RAD21 IAA vs RAD21 control | 84 | 21,620 | 77 | 0.36% |
| WAPL IAA vs WAPL control | 0 | 21,620 | 0 | 0% |

### CTCF

CTCF depletion produced a strong RNA-level response: 1,981 significant DGE genes and 6,795 significant APA sites. This is consistent with a broad transcriptional and 3′-end response. However, 31.4% of tested PAS being significant is extensive and should be verified using:

- replicate clustering and pair behavior in the C4 PCA;
- effect-size and adjusted-p-value distributions;
- distal versus proximal shift counts;
- internal-priming and PAS-localization summaries;
- representative loci in IGV.

A real global response and a systematic compositional or batch difference can both produce many significant PAS; the current counts alone cannot distinguish them.

### RAD21

RAD21 depletion produced a modest but detectable response: 84 DGE genes and 77 significant APA sites. This resembles a relatively limited or specific RNA-level response rather than a technical failure.

### WAPL

No significant DGE genes or APA sites were detected for WAPL IAA versus control. This does not establish absence of a biological WAPL effect. Possible explanations include:

- genuinely weak RNA-level effects at the tested time point;
- weak or failed WAPL degradation;
- insufficient treatment duration;
- reduced statistical power caused by variable IAA libraries;
- one outlying replicate;
- sample-label, pairing or batch problems.

`WAPL_IAA_R3` has only 14.37 million STAR-input reads and 11.37 million unique C0 reads, approximately half the depth of the other WAPL IAA libraries. The WAPL IAA group also has the greatest depth variability. Its PCA position, sample correlations, size factor and orthogonal WAPL-degradation evidence should be reviewed before interpreting zero discoveries as a biological negative.

## All configured contrasts

| Contrast | Design | Significant DGE genes | Significant APA-A sites |
|---|---|---:|---:|
| CTCF IAA vs CTCF control | paired | 1,981 | 6,795 |
| RAD21 control vs CTCF control | unpaired | 4,798 | 5,128 |
| RAD21 IAA vs CTCF control | unpaired | 5,006 | 3,899 |
| WAPL control vs CTCF control | unpaired | 2,251 | 1,501 |
| WAPL IAA vs CTCF control | unpaired | 1,079 | 234 |
| RAD21 control vs CTCF IAA | unpaired | 4,011 | 8,670 |
| RAD21 IAA vs CTCF IAA | unpaired | 3,829 | 6,809 |
| WAPL control vs CTCF IAA | unpaired | 1,420 | 4,454 |
| WAPL IAA vs CTCF IAA | unpaired | 569 | 2,848 |
| RAD21 IAA vs RAD21 control | paired | 84 | 77 |
| WAPL control vs RAD21 control | unpaired | 1,589 | 963 |
| WAPL IAA vs RAD21 control | unpaired | 74 | 10 |
| WAPL control vs RAD21 IAA | unpaired | 1,886 | 1,233 |
| WAPL IAA vs RAD21 IAA | unpaired | 75 | 42 |
| WAPL IAA vs WAPL control | paired | 0 | 0 |

The 12 cross-target comparisons involve disjoint subjects and different degron backgrounds. They can contain target, clone, baseline-state and batch differences and should be treated as exploratory. They are not direct estimates of the IAA effect. The large cross-target counts therefore do not supersede the three paired within-target results.

## Samples and patterns requiring attention

1. **WAPL IAA R3:** lowest depth in the experiment; review MultiQC, PCA, normalized counts and sample correlations.
2. **WAPL IAA group:** highest depth variability; zero paired discoveries may reflect limited power.
3. **RAD21 IAA R1/R2:** lowest unique-mapping fractions and relatively high multimapping; review repetitive/rRNA/low-complexity indicators.
4. **CTCF APA response:** very large significant fraction; validate replicate agreement and directional/genomic distribution.
5. **Cross-target contrasts:** do not interpret as degron treatment effects without a model explicitly addressing target background and batch.

## Reports and files worth reviewing

Alpha.9 creates a final HTML file, but it is mainly a module/run-status report rather than a comprehensive scientific report. It does not integrate STAR QC, DGE direction, APA shifts, PCPA, enrichment and browser evidence into one interpretation. Alpha.10 is intended to provide that integrated report.

Use the following files from the completed alpha.9 run, in roughly this order:

1. `01_qc/multiqc/multiqc_report.html` — read quality, trimming, mapping and sample outliers.
2. `02_alignment/protocol_orientation.tsv` — QuantSeq REV strand compatibility.
3. `05_gene_expression/GRCm39/C4_primary_deseq2/C4_vst_pca.pdf` — replicate clustering and outlier detection.
4. `05_gene_expression/GRCm39/C4_primary_deseq2/result_index.tsv` — all DGE contrasts.
5. `05_gene_expression/GRCm39/C4_primary_deseq2/<contrast>.deseq2.tsv` and `<contrast>.MA.pdf` — effect sizes and significance.
6. `05_gene_expression/GRCm39/C4_vs_C5_library_correlations.tsv` and `C4_vs_C5_large_discrepancies.tsv` — comparison of 3′-end gene sums with conventional exon counting.
7. `06_apa_a_mcell2019/GRCm39/dexseq/result_index.tsv` — all APA-A contrasts.
8. `06_apa_a_mcell2019/GRCm39/dexseq/<contrast>.dexseq.tsv` — PAS-level APA statistics.
9. `06_apa_a_mcell2019/GRCm39/dexseq/<contrast>.apa_shift.tsv` — gene-level distal/proximal/no-shift classifications.
10. `06_apa_a_mcell2019/GRCm39/candidate_pcpa.tsv` — significant differential candidate PCPA events.
11. `04_active_pas/GRCm39/active_pas_catalog.tsv` and `count_universe_audit.json` — PAS annotation and count-universe validation.
12. `10_reports/IGV_session.xml`, `UCSC_trackDb.txt` and `09_tracks/track_normalization.tsv` — genome-browser review and normalization provenance.

In the DGE index, `significant` means adjusted p-value below the configured FDR and does not impose an additional fold-change threshold. In the APA index, `significant_sites` means DEXSeq-adjusted p-value below the FDR. Differential PCPA output additionally requires the configured minimum absolute `delta_PAU`.

## Remaining high-priority checks

Before a firm biological conclusion, perform the following:

1. inspect the C4 VST PCA, especially `WAPL_IAA_R3` and the control/IAA pairs;
2. review MultiQC for the low-depth and lower-mapping samples;
3. inspect C4 size factors and C4/C5 correlations;
4. extract the per-sample C0/C1/C1S/C2/C2R funnel;
5. count distal, proximal, no-shift and not-classifiable genes for the three paired contrasts;
6. summarize terminal-exon, intronic, internal-exonic, downstream, antisense, intergenic and ambiguous PAS classes;
7. inspect representative strong DGE, APA and candidate PCPA loci in IGV;
8. confirm CTCF, RAD21 and WAPL degradation independently, for example by protein-level assay.

Corrected command for the per-sample end-processing funnel:

```bash
RESULTS="$HOME/Analysis/mesc_degron_run8_rna_ends2tracks/results_alpha9_pas_v1"
PYTHON="/opt/conda_envs/rna_ends2tracks-0.1.0a9.post2/bin/python"

"$PYTHON" -c 'import glob,json,os,sys
print("sample\tC0\tC1\tC1S\tC2\tC2R\tC2_retained_pct")
for p in sorted(glob.glob(sys.argv[1]+"/03_exact_ends/GRCm39/*/end_audit.json")):
 d=json.load(open(p))
 c1=int(d.get("C1",0)); c2=int(d.get("C2",0))
 retained=100*c2/c1 if c1 else 0
 print("{}\t{}\t{}\t{}\t{}\t{}\t{:.2f}".format(os.path.basename(os.path.dirname(p)),d.get("C0",0),c1,d.get("C1S",0),c2,d.get("C2R",0),retained))' \
"$RESULTS"
```

The earlier diagnostic `SyntaxError` came from an invalid escaped expression in the display command. It did not affect the completed workflow or any result file.

## Preliminary conclusion

The sequencing experiment overall worked technically: it has deep coverage, good unique mapping, strong and consistent QuantSeq REV strand specificity, a large filtered exact-end universe, high C2-to-C3 assignment, valid PAS-count invariants and complete DGE/APA outputs.

The preliminary biological picture is a strong CTCF-associated transcriptional and APA response, a modest RAD21 response, and no detected WAPL response. The CTCF result requires validation because of its unusually broad APA signal, while the WAPL result requires particular caution because variable and lower-depth IAA libraries may have reduced statistical power. PCA, per-sample C2 retention, PAS-location composition and orthogonal degron-efficiency evidence are the most important next checks.
