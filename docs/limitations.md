# Limitations and interpretation

- Lexogen QuantSeq REV V1/V2 single-end and paired-end no-UMI profiles are accepted. A project uses one layout consistently; mixed SE/PE projects are rejected. In PE mode both mates are aligned and contribute their aligned blocks to conventional coverage (normalized per mapped pair), while only R1 defines the cleavage coordinate. The all-read track is not an inferred insert-span track. UMI protocols remain unsupported.
- PCR duplicates cannot be identified reliably without UMIs. Duplicate flags are retained, so amplification bias may remain.
- Exact cleavage coordinates depend on alignment/CIGAR interpretation. End-defining clipped reads are excluded to avoid false precision.
- Internal-priming masking is sequence-rule-based. Atlas rescue improves sensitivity but may rescue false sites; `core_plus_rescue` is a sensitivity mode.
- Active PAS are project-specific. Adding samples can change pooled discovery and therefore C3/C4; comparisons across independently discovered PAS universes require harmonization.
- C4 DGE assumes most genes are not changing in one direction. DESeq2/robust tracks are relative cohort scaling, not absolute calibration.
- C4 and conventional C5 quantify different signal definitions. Large discrepancies require biological/annotation review rather than automatic substitution.
- DEXSeq identifies differential PAS usage, not necessarily a change in total gene output.
- Intragenic PAS are candidate PCPA/premature termination events. Confirmation needs orthogonal evidence such as nascent transcription, long-read RNA or validated isoform assays.
- Ambiguous multi-gene PAS are excluded from statistical matrices by default. This is conservative but can reduce sensitivity in overlapping loci.
- APA-B uses a separately pinned PolyAseqTrap/DeepIP environment and is validation-scope gated. The audited `biolserv` GRCh38/GRCm39 QuantSeq REV V2 single-end scope is accepted; other assemblies or protocols, including paired-end, require their own real-data acceptance. Its results are independent and compared only after both branches complete; it is not a fallback for APA-A.
- APA-A2 corrects effect estimation and shift classification but intentionally shares APA-A's condition-blind C3 site catalog. It is an independent statistical run, not an independent PAS-discovery method. Retain APA-B when discovery-method independence is required.
- Legacy accepted APA-B manifests created before the protocol-list field are interpreted narrowly as QuantSeq REV V2 SE only. They do not authorize PE execution; each PE protocol requires a reviewed real-PE canary.
- Human GRCh38/GENCODE v42 and mouse GRCm39/GENCODE vM31 are the initial validated profiles. Other assemblies/releases require a new synchronized reference and PAS-atlas profile.
