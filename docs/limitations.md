# Limitations and interpretation

- Only Lexogen QuantSeq REV V1/V2 single-end, no-UMI profiles are accepted. Paired-end and UMI protocols need separate validated coordinate/count policies.
- PCR duplicates cannot be identified reliably without UMIs. Duplicate flags are retained, so amplification bias may remain.
- Exact cleavage coordinates depend on alignment/CIGAR interpretation. End-defining clipped reads are excluded to avoid false precision.
- Internal-priming masking is sequence-rule-based. Atlas rescue improves sensitivity but may rescue false sites; `core_plus_rescue` is a sensitivity mode.
- Active PAS are project-specific. Adding samples can change pooled discovery and therefore C3/C4; comparisons across independently discovered PAS universes require harmonization.
- C4 DGE assumes most genes are not changing in one direction. DESeq2/robust tracks are relative cohort scaling, not absolute calibration.
- C4 and conventional C5 quantify different signal definitions. Large discrepancies require biological/annotation review rather than automatic substitution.
- DEXSeq identifies differential PAS usage, not necessarily a change in total gene output.
- Intragenic PAS are candidate PCPA/premature termination events. Confirmation needs orthogonal evidence such as nascent transcription, long-read RNA or validated isoform assays.
- Ambiguous multi-gene PAS are excluded from statistical matrices by default. This is conservative but can reduce sensitivity in overlapping loci.
- APA-B is an external pilot-gated method. Its results are independent and compared only after both branches complete; it is not a fallback for APA-A.
- Human GRCh38/GENCODE v42 and mouse GRCm39/GENCODE vM31 are the initial validated profiles. Other assemblies/releases require a new synchronized reference and PAS-atlas profile.
