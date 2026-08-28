# 11. APA-B and comparison with APA-A

APA-B answers the same biological question with a deliberately different discovery method. It clusters genome-wide QuantSeq REV 3′ endpoints with PolyAseqTrap, filters sequence-based internal-priming candidates with DeepIP, and tests PAS usage using DRIMSeq/stageR. APA-A remains the transparent Mcell2019-style primary analysis.

Both methods can detect PAS inside gene bodies. An internal-exonic or intronic PAS upstream of the normal terminal PAS is reported as a candidate premature cleavage/polyadenylation (PCPA) event. The assay does not directly measure polymerase occupancy, so the report must not call it proven premature transcription termination without independent evidence.

APA-B is disabled until the shared installation passes a synthetic and real-data pilot. This is a one-time site validation, not something every project user repeats. After acceptance, the normal project still has only two inputs: `config.conf` and `samplesheet.csv`.

The two PAS catalogs are never combined. `08_apa_comparison/` matches nearby sites only after independent statistics and reports:

- sites detected by both methods;
- APA-A-only and APA-B-only sites;
- concordant and discordant directions of PAS-usage change;
- proximity between method-specific coordinates.

Agreement is supportive evidence, not a replacement for inspecting read coverage, replicate consistency, gene structure, and internal-priming evidence. See [the installation and pilot contract](POLYASEQTRAP_ADAPTER_CONTRACT.md) for exact pins, commands, outputs, and limitations.
