# Alpha.3 deployment-readiness review

> Historical alpha.3 validation record; it is not the alpha.6 run or installation guide.

**Audience:** workflow owner and `biolserv` administrator  
**Review date:** 2026-08-25  
**Scope:** no-UMI Lexogen QuantSeq REV single-end processing, human GRCh38 and mouse GRCm39, the shared Conda installation, and reuse of the audited STAR/reference assets.

## Direct answer

The isolated alpha.3 environment is ready for a versioned **server canary installation**, but the workflow is not yet ready to be described as production validated. The package specification solved on `biolserv` with Mamba 2.8.0 (386 packages; approximately 450 MB). Final pre-tag commit `b1133268ea852d528e25b9910f3f80d708743154` passed both [CI](https://github.com/MichalGd/rna_ends2tracks/actions/runs/32838086483) and [Complete environment](https://github.com/MichalGd/rna_ends2tracks/actions/runs/32838086350).

The research review found and corrected one pre-tag issue: BBDuk was invoked with `qtrim=t`, which BBTools interprets as trimming both ends. Because the APA-defining QuantSeq REV coordinate is the 5-prime end of Read 1, alpha.3 now uses `qtrim=r`, preserving that coordinate. A unit assertion and a real BBDuk CI smoke test protect this command.

Proceed in this order: tag alpha.3, create a new environment prefix, install the tagged wheel or immutable tagged source, validate the two reference manifests, and run synthetic/small real-data canaries. Do not enable APA-B or call the deployment production-ready yet.

## Evidence and decisions

### Environment resolution and isolation

- Bioconda currently requires `conda-forge` at higher priority than `bioconda`, strict channel priority, and no `defaults`; `defaults` was removed from its recommendation in August 2024 ([Bioconda usage guidance](https://bioconda.github.io/)). The repository now encodes `conda-forge`, `bioconda`, and `nodefaults`.
- Mamba accepts configuration through environment variables as well as RC files and command-line options ([Mamba configuration](https://mamba.readthedocs.io/en/latest/user_guide/configuration.html)). `MAMBA_CHANNEL_PRIORITY=strict` with `--no-rc` is therefore appropriate for Mamba 2.8.0 and does not alter any user's `.condarc`.
- The official setup action supports `environment-file`, environment caching, and `micromamba-shell` ([setup-micromamba documentation](https://github.com/mamba-org/setup-micromamba)). The complete Linux CI created the full environment and verified the required Python, R, and command-line packages.
- The new environment must remain a separate, versioned prefix. Existing `rnaseq2tracks` and `cutnrun2tracks` environments can contribute downloaded package-cache reuse, but their installed binaries must not be combined at runtime with this workflow.

### Lexogen REV raw-read processing

- Lexogen's January 2026 REV V2 guide states that Read 1 represents cDNA, is opposite to the genomic-reference orientation, and pinpoints the exact 3-prime end. It recommends FastQC, removal of adapters, poly(A)/poly(T), low-quality sequence, and reads shorter than 20 nt, followed by splice-aware STAR alignment ([QuantSeq REV V2 user guide](https://www.lexogen.com/wp-content/uploads/2026/01/225UG675V0102_QuantSeq-REV-V2-3-mRNA-Seq_2026-01-19.pdf)).
- Lexogen's current FAQ says its REV V2 pipeline uses single-read data and asks users to upload only Read 1 ([REV V2 analysis FAQ](https://faqs.lexogen.com/faq/how-can-i-analyze-my-quantseq-rev-v2-data)). This supports the workflow's explicit SE/R1-only contract.
- The older published Lexogen pipeline gives the concrete BBDuk settings `ktrim=r`, `qtrim=r`, `trimq=10`, and `minlength=20`, and uses STAR followed by reverse-stranded counting ([Lexogen integrated pipeline guide](https://www.lexogen.com/wp-content/uploads/2019/12/015UG108V0200-QuantSeq-Data-Analysis-Pipeline_2018-10-18.pdf)).
- BBTools documents `qtrim=r` as right-end trimming and `qtrim=rl` as both-end trimming ([BBDuk guide](https://github.com/BioInfoTools/BBMap/blob/master/docs/guides/BBDukGuide.txt)). The source parser also treats boolean true as enabling both left and right trimming ([BBTools Parser source](https://raw.githubusercontent.com/bbushnell/BBTools/master/current/parse/Parser.java)). For exact-end APA, the right-only choice is therefore consequential rather than cosmetic.
- No UMI extraction or deduplication is performed. The `umi: false` protocol profile and retention of coordinate-identical primary reads are appropriate for the requested non-UMI libraries. The kit's unique dual sample indices are demultiplexing indices, not molecular UMIs.

### Existing STAR indexes

- STAR defines `sjdbOverhang` as ideally read length minus one ([STAR parameter definitions](https://github.com/alexdobin/STAR/blob/master/source/parametersDefault)). The STAR author has clarified that indexes built for longer reads work for shorter reads and that one default/long-read index is normally sufficient for reads over 50 nt; very short reads such as 25 nt deserve a separate index ([STAR support discussion](https://groups.google.com/g/rna-star/c/x60p1C-pGbc)).
- Consequently, the audited `sjdbOverhang=149` indexes are reasonable canary candidates for 50, 75, 100, and 150 nt libraries. Separate indexes for every read length are not justified. Retain the workflow warning when the declared maximum read length differs, and compare a representative 50 nt sample only if that length is actually used.
- Structural compatibility has passed: STAR-required files exist and contig names/lengths match the selected FASTA indexes. Runtime loading by STAR 2.7.11b and a small alignment remain mandatory because structural checks alone do not prove runtime compatibility or biological equivalence.

## Remaining gates

| Gate | Status | Required action |
|---|---|---|
| Alpha.3 environment solve | Passed | Preserve the solved specification and create only a new prefix. |
| Full Linux package integration | Passed on final pre-tag commit | Preserve the successful Actions run with the release record. |
| REV coordinate-preserving trimming | Passed | Unit tests and the real BBDuk environment smoke test passed on `b1133268`. |
| Human and mouse reference manifests | Audited assets available | Create release manifests with exact absolute paths and SHA-256 provenance. |
| Chromosome sizes | Full files generated in audit area | Publish immutable copies alongside each release manifest; do not use the subset-like legacy files. |
| PAS atlas | Missing | APA-A can run without known-PAS annotation, but interpretation is reduced; obtain assembly-matched resources when feasible. |
| APA-A scientific validation | Pending | Run synthetic exact-end/internal-priming/PCPA truth sets and a small REV canary. |
| APA-B PolyAseqTrap/DeepIP | Intentionally blocked | Install separately, pin code/models, and complete the documented REV/no-UMI pilot before setting `pilot_accepted: true`. |
| Shared immutability | Pending | After validation, make the versioned software prefix and manifests administrator-owned and non-writable to ordinary users. Keep project output paths writable. |

## Installation decision

Creating the new Conda prefix is safe with respect to existing analyses because it does not modify their environments or references. It can consume package-cache I/O, CPU, and network bandwidth, but the server audit showed ample capacity and no active analysis jobs. Never update an existing shared environment in place, replace reference files in place, or publish an unversioned `current` launcher during this canary.

The first installation should expose DGE and APA-A for validation while leaving APA-B disabled. A successful package installation is not sufficient: the release is accepted only after command/version checks, both reference validations, a dry run, a tiny human canary, a tiny mouse canary, and inspection of orientation and exact-end outputs.

## Material limitations

Lexogen's current guidance describes general preprocessing and exact REV end behavior but does not validate this repository's complete APA statistics. The existing STAR indexes have been structurally audited but are still owner-writable and have not yet been loaded by the new runtime. PolyAseqTrap/DeepIP and an assembly-matched PAS atlas are not installed. Candidate intragenic cleavage sites can be described as candidate PCPA events consistent with premature termination, but QuantSeq alone cannot prove downstream RNA-polymerase termination.
