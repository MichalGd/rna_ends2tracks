from __future__ import annotations

import csv
import gzip
import json
import re
from collections import defaultdict
from collections.abc import Iterable
from dataclasses import dataclass, field
from pathlib import Path

import pysam

from .config import RunPlan, signature_for
from .external import event, require_tools, run
from .receipts import receipt_valid, write_receipt

PAS_MOTIFS = ("AATAAA", "ATTAAA", "TATAAA", "AGTAAA", "AAGAAA", "AATATA", "AATACA", "CATAAA")


@dataclass
class Gene:
    gene_id: str
    chrom: str
    start: int
    end: int
    strand: str
    exons: list[tuple[int, int]] = field(default_factory=list)
    cds: list[tuple[int, int]] = field(default_factory=list)
    utr: list[tuple[int, int]] = field(default_factory=list)


def reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1]


def transcript_end(read: pysam.AlignedSegment) -> tuple[int, str, bool]:
    """Return zero-based cleavage base, transcript strand, and end-defining soft-clip flag."""
    if read.is_reverse:
        clipped = bool(read.cigartuples and read.cigartuples[-1][0] == 4)
        return read.reference_end - 1, "+", clipped
    clipped = bool(read.cigartuples and read.cigartuples[0][0] == 4)
    return read.reference_start, "-", clipped


def _attributes(text: str) -> dict[str, str]:
    result: dict[str, str] = {}
    for item in text.strip().strip(";").split(";"):
        if not item.strip():
            continue
        key, _, value = item.strip().partition(" ")
        result[key] = value.strip().strip('"')
    return result


def load_genes(gtf_path: str) -> tuple[dict[str, Gene], dict[tuple[str, str, int], set[str]]]:
    opener = gzip.open if gtf_path.endswith(".gz") else open
    genes: dict[str, Gene] = {}
    with opener(gtf_path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] not in {"gene", "exon", "CDS", "three_prime_utr", "UTR"}:
                continue
            attrs = _attributes(fields[8])
            gene_id = attrs.get("gene_id")
            if not gene_id:
                continue
            chrom, feature, start, end, strand = fields[0], fields[2], int(fields[3]) - 1, int(fields[4]), fields[6]
            if gene_id not in genes:
                genes[gene_id] = Gene(gene_id, chrom, start, end, strand)
            gene = genes[gene_id]
            gene.start = min(gene.start, start)
            gene.end = max(gene.end, end)
            if feature == "exon":
                gene.exons.append((start, end))
            elif feature == "CDS":
                gene.cds.append((start, end))
            elif feature in {"three_prime_utr", "UTR"}:
                gene.utr.append((start, end))
    bins: dict[tuple[str, str, int], set[str]] = defaultdict(set)
    bin_size = 100_000
    for gene in genes.values():
        for bin_number in range(gene.start // bin_size, max(gene.start, gene.end - 1) // bin_size + 1):
            bins[(gene.chrom, gene.strand, bin_number)].add(gene.gene_id)
    return genes, bins


def load_known_pas(path: str | None) -> dict[tuple[str, str], list[int]]:
    result: dict[tuple[str, str], list[int]] = defaultdict(list)
    if not path:
        return result
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if not line or line.startswith(("#", "track", "browser")):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                continue
            strand = fields[5] if len(fields) > 5 and fields[5] in {"+", "-"} else "."
            result[(fields[0], strand)].append(int(fields[1]))
    for values in result.values():
        values.sort()
    return result


def nearest_known(chrom: str, strand: str, position: int, known: dict[tuple[str, str], list[int]]) -> int | None:
    import bisect
    values = known.get((chrom, strand), []) + known.get((chrom, "."), [])
    if not values:
        return None
    values.sort()
    index = bisect.bisect_left(values, position)
    candidates = values[max(0, index - 1): index + 1]
    return min(abs(position - value) for value in candidates)


def _overlaps(intervals: Iterable[tuple[int, int]], position: int) -> bool:
    return any(start <= position < end for start, end in intervals)


def annotate_site(chrom: str, strand: str, position: int, genes: dict[str, Gene], bins: dict[tuple[str, str, int], set[str]]) -> tuple[str, str]:
    compatible: list[tuple[str, str]] = []
    for gene_id in bins.get((chrom, strand, position // 100_000), set()):
        gene = genes[gene_id]
        if not gene.start <= position < gene.end:
            continue
        terminal_exons = [interval for interval in gene.exons if interval[1] == gene.end] if strand == "+" else [interval for interval in gene.exons if interval[0] == gene.start]
        if _overlaps(gene.utr, position):
            if gene.cds:
                is_three_prime = position >= max(end for _, end in gene.cds) if strand == "+" else position < min(start for start, _ in gene.cds)
            else:
                is_three_prime = _overlaps(terminal_exons, position)
            feature = "terminal_3UTR" if is_three_prime else "internal_UTR"
        elif _overlaps(gene.cds, position):
            feature = "internal_exon_CDS"
        elif _overlaps(gene.exons, position):
            feature = "terminal_exon" if _overlaps(terminal_exons, position) else "internal_exon"
        else:
            feature = "intron"
        compatible.append((gene_id, feature))
    if len(compatible) == 1:
        return compatible[0]
    if len(compatible) > 1:
        return "", "ambiguous_multi_gene"
    return "", "intergenic"


def _sequence_windows(fasta: pysam.FastaFile, chrom: str, position: int, strand: str, downstream: int = 20) -> tuple[str, str]:
    chrom_length = fasta.get_reference_length(chrom)
    if strand == "+":
        down = fasta.fetch(chrom, min(chrom_length, position + 1), min(chrom_length, position + 1 + downstream))
        up = fasta.fetch(chrom, max(0, position - 50), max(0, position - 10))
    else:
        down = reverse_complement(fasta.fetch(chrom, max(0, position - downstream), max(0, position)))
        up = reverse_complement(fasta.fetch(chrom, min(chrom_length, position + 10), min(chrom_length, position + 50)))
    return down.upper(), up.upper()


def _longest_a(sequence: str) -> int:
    return max((len(match) for match in re.findall(r"A+", sequence)), default=0)


def apa_a(plan: RunPlan, results: Path, script_root: Path, dry_run: bool = False, force: bool = False) -> None:
    module_dir = results / "04_apa_a_repository"
    log_dir = results / "provenance" / "logs"
    module_dir.mkdir(parents=True, exist_ok=True)
    bams = [results / "02_alignment" / s["sample_id"] / f"{s['sample_id']}.bam" for s in plan.samples]
    if dry_run:
        event(log_dir, "apa_a", "dry_run", "Would extract exact REV ends, cluster sites, annotate, count, and run DEXSeq")
        return
    missing = [str(path) for path in bams if not path.is_file()]
    if missing:
        raise RuntimeError("APA-A requires completed sample BAMs: " + ", ".join(missing))
    signature = signature_for([*bams, plan.reference["fasta"], plan.reference["gtf"]], {
        "module": "apa_a", "parameters": plan.project.get("apa_a", {}), "design": plan.project["design"],
        "samples": plan.samples, "contrasts": plan.contrasts, "reporting": plan.project.get("reporting", {}),
    })
    catalog_path = module_dir / "pas_catalog.tsv"
    counts_path = module_dir / "pas_counts.tsv"
    exact_path = module_dir / "exact_ends.tsv.gz"
    pcpa_catalog_path = module_dir / "pcpa_candidate_catalog.tsv"
    pcpa_path = module_dir / "candidate_pcpa.tsv"
    extraction_audit_path = module_dir / "end_extraction_audit.json"
    stats_index = module_dir / "dexseq" / "result_index.tsv"
    if not force and receipt_valid(module_dir, signature):
        event(log_dir, "apa_a", "skipped", "Valid matching receipt")
        return

    settings = plan.project.get("apa_a", {})
    min_reads = int(settings.get("min_reads", 5))
    min_samples = int(settings.get("min_samples", 2))
    cluster_gap = int(settings.get("cluster_gap_nt", 24))
    known_distance = int(settings.get("max_known_pas_distance_nt", 24))
    counts: dict[tuple[str, str, int], dict[str, int]] = defaultdict(lambda: defaultdict(int))
    audit = {"mapped_records": 0, "unique_primary_records": 0, "end_soft_clipped_records": 0,
             "duplicate_flagged_unique_primary_records": 0, "primary_counted_records": 0}
    exact_path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(exact_path, "wt", encoding="utf-8", newline="") as out:
        writer = csv.writer(out, delimiter="\t", lineterminator="\n")
        writer.writerow(["chrom", "start", "end", "strand", "sample_id", "read_name", "mapping_class", "duplicate_flag", "nh", "end_soft_clipped", "primary_eligible"])
        for sample, bam_path in zip(plan.samples, bams):
            with pysam.AlignmentFile(bam_path, "rb") as bam:
                for read in bam.fetch(until_eof=True):
                    if plan.project.get("protocol", {}).get("library_layout", "SE") == "PE" and not read.is_read1:
                        audit["non_end_defining_mate_records"] += 1
                        continue
                    if read.is_unmapped:
                        continue
                    audit["mapped_records"] += 1
                    nh = read.get_tag("NH") if read.has_tag("NH") else 0
                    position, strand, clipped = transcript_end(read)
                    mapping_class = "secondary" if read.is_secondary else "supplementary" if read.is_supplementary else "primary_unique" if nh == 1 else "primary_multimapping"
                    eligible = not read.is_secondary and not read.is_supplementary and nh == 1
                    writer.writerow([read.reference_name, position, position + 1, strand, sample["sample_id"], read.query_name,
                                     mapping_class, int(read.is_duplicate), nh, int(clipped), int(eligible and not clipped)])
                    if eligible:
                        audit["unique_primary_records"] += 1
                        if read.is_duplicate:
                            audit["duplicate_flagged_unique_primary_records"] += 1
                    if clipped:
                        audit["end_soft_clipped_records"] += 1
                    if eligible and not clipped:
                        counts[(read.reference_name, strand, position)][sample["sample_id"]] += 1
                        audit["primary_counted_records"] += 1
    extraction_audit_path.write_text(json.dumps(audit, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    clusters: list[list[tuple[str, str, int]]] = []
    for (chrom, strand), keys in _group_positions(counts).items():
        current: list[tuple[str, str, int]] = []
        previous: int | None = None
        for key in sorted(keys, key=lambda item: item[2]):
            if previous is None or key[2] - previous <= cluster_gap:
                current.append(key)
            else:
                clusters.append(current)
                current = [key]
            previous = key[2]
        if current:
            clusters.append(current)

    genes, bins = load_genes(plan.reference["gtf"])
    known = load_known_pas(plan.reference.get("pas_atlas"))
    fasta = pysam.FastaFile(plan.reference["fasta"])
    sample_ids = [sample["sample_id"] for sample in plan.samples]
    catalog: list[dict[str, object]] = []
    matrix: list[dict[str, object]] = []
    for members in clusters:
        member_count = {sample: sum(counts[key].get(sample, 0) for key in members) for sample in sample_ids}
        total = sum(member_count.values())
        supporting = sum(value > 0 for value in member_count.values())
        if total < min_reads or supporting < min_samples:
            continue
        summit = min(members, key=lambda key: (-sum(counts[key].values()), key[2]))
        chrom, strand, position = summit
        down, up = _sequence_windows(fasta, chrom, position, strand)
        a_fraction = down.count("A") / len(down) if down else 0.0
        longest_a = _longest_a(down)
        motif = next((motif for motif in PAS_MOTIFS if motif in up), "")
        distance = nearest_known(chrom, strand, position, known)
        a_rich = a_fraction >= float(settings.get("internal_priming_a_fraction", 0.70)) or longest_a >= int(settings.get("internal_priming_a_run", 6))
        rescued = a_rich and (bool(motif) or (distance is not None and distance <= known_distance))
        confidence = "rescued_a_rich" if rescued else "probable_internal_priming" if a_rich else "high_confidence"
        gene_id, feature = annotate_site(chrom, strand, position, genes, bins)
        pas_id = f"APA_A_{chrom}_{position + 1}_{'P' if strand == '+' else 'M'}"
        catalog.append({
            "pas_id": pas_id, "chrom": chrom, "start": position, "end": position + 1, "strand": strand,
            "summit_1based": position + 1, "cluster_start": min(x[2] for x in members),
            "cluster_end": max(x[2] for x in members) + 1, "member_positions": ",".join(str(x[2] + 1) for x in members),
            "total_reads": total, "supporting_samples": supporting, "downstream_a_fraction": round(a_fraction, 4),
            "longest_a_run": longest_a, "upstream_pas_motif": motif, "known_pas_distance": "" if distance is None else distance,
            "confidence": confidence, "gene_id": gene_id, "feature_class": feature,
        })
        matrix.append({"pas_id": pas_id, **member_count})
    fasta.close()
    _write_rows(catalog_path, catalog)
    _write_rows(counts_path, matrix)
    _write_pcpa(pcpa_catalog_path, catalog, genes)

    require_tools(["Rscript"])
    run([
        "Rscript", str(script_root / "R" / "dexseq_all_pairs.R"), "--counts", str(counts_path),
        "--catalog", str(catalog_path), "--samples", str(results / "00_metadata" / "validated_samples.tsv"),
        "--contrasts", str(results / "00_metadata" / "contrasts.tsv"), "--outdir", str(module_dir / "dexseq"),
        "--min-count", str(min_reads), "--design", str(plan.project["design"]),
    ], log_dir / "apa_a" / "dexseq.log", False)
    _filter_significant_pcpa(
        pcpa_catalog_path, catalog_path, stats_index, pcpa_path,
        float(plan.project.get("reporting", {}).get("fdr", 0.05)),
        float(plan.project.get("reporting", {}).get("min_abs_delta_pau", 0.10)),
    )
    write_receipt("apa_a", module_dir, signature, [catalog_path, counts_path, exact_path, extraction_audit_path, pcpa_catalog_path, pcpa_path, stats_index], ["rna-ends2tracks", "apa-a"])
    event(log_dir, "apa_a", "completed", f"Discovered {len(catalog)} PAS clusters")


def _group_positions(counts: dict[tuple[str, str, int], dict[str, int]]) -> dict[tuple[str, str], list[tuple[str, str, int]]]:
    result: dict[tuple[str, str], list[tuple[str, str, int]]] = defaultdict(list)
    for key in counts:
        result[(key[0], key[1])].append(key)
    return result


def _write_rows(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise RuntimeError(f"No rows passed filters for {path.name}")
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def _write_pcpa(path: Path, catalog: list[dict[str, object]], genes: dict[str, Gene]) -> None:
    terminal_genes = {str(row["gene_id"]) for row in catalog
                      if str(row["feature_class"]).startswith("terminal") and row["gene_id"]
                      and row["confidence"] in {"high_confidence", "rescued_a_rich"}}
    rows: list[dict[str, object]] = []
    eligible_features = {"intron", "internal_exon", "internal_exon_CDS"}
    for row in catalog:
        gene_id = str(row["gene_id"])
        if gene_id not in terminal_genes or row["feature_class"] not in eligible_features:
            continue
        if row["confidence"] not in {"high_confidence", "rescued_a_rich"}:
            continue
        gene = genes[gene_id]
        position = int(row["start"])
        terminal = gene.end - 1 if gene.strand == "+" else gene.start
        upstream = position < terminal if gene.strand == "+" else position > terminal
        if not upstream:
            continue
        consequence = "coding_truncating_intronic_PCPA" if row["feature_class"] in {"intron", "internal_exon_CDS"} else "upstream_exonic_termination"
        rows.append({
            "pas_id": row["pas_id"], "gene_id": gene_id, "chrom": row["chrom"], "start": position,
            "end": row["end"], "strand": row["strand"], "feature_class": row["feature_class"],
            "consequence": consequence, "confidence": row["confidence"],
            "interpretation": "candidate PCPA consistent with premature transcription termination",
        })
    headers = ["pas_id", "gene_id", "chrom", "start", "end", "strand", "feature_class", "consequence", "confidence", "interpretation"]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=headers, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def _filter_significant_pcpa(candidates_path: Path, catalog_path: Path, index_path: Path, output_path: Path, fdr: float, min_delta: float) -> None:
    with candidates_path.open(encoding="utf-8", newline="") as handle:
        candidates = {row["pas_id"]: row for row in csv.DictReader(handle, delimiter="\t")}
    with catalog_path.open(encoding="utf-8", newline="") as handle:
        catalog = {row["pas_id"]: row for row in csv.DictReader(handle, delimiter="\t")}
    output_rows: list[dict[str, object]] = []
    with index_path.open(encoding="utf-8", newline="") as handle:
        indexes = list(csv.DictReader(handle, delimiter="\t"))
    for index in indexes:
        with Path(index["result_file"]).open(encoding="utf-8", newline="") as handle:
            results = list(csv.DictReader(handle, delimiter="\t"))
            tested_terminal_genes = {catalog[row["pas_id"]]["gene_id"] for row in results
                                     if row["pas_id"] in catalog and catalog[row["pas_id"]]["feature_class"].startswith("terminal")}
            for result in results:
                if result["pas_id"] not in candidates or result.get("padj", "NA") in {"", "NA"}:
                    continue
                if candidates[result["pas_id"]]["gene_id"] not in tested_terminal_genes:
                    continue
                if float(result["padj"]) <= fdr and abs(float(result.get("delta_PAU", 0))) >= min_delta:
                    output_rows.append({
                        **candidates[result["pas_id"]], "contrast_id": index["contrast_id"],
                        "padj": result["padj"], "delta_PAU": result.get("delta_PAU", ""),
                    })
    base_headers = ["pas_id", "gene_id", "chrom", "start", "end", "strand", "feature_class", "consequence", "confidence", "interpretation"]
    headers = [*base_headers, "contrast_id", "padj", "delta_PAU"]
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=headers, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(output_rows)
