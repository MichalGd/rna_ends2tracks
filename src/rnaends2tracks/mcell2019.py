from __future__ import annotations

import bisect
import csv
import gzip
from collections import defaultdict
from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from .apa_a import Gene, load_genes


@dataclass(frozen=True)
class Site:
    chrom: str
    strand: str
    summit: int
    start: int
    end: int


def defining_end_is_clipped(cigartuples: list[tuple[int, int]] | None, is_reverse: bool) -> bool:
    """CIGAR is reference-left to reference-right; REV cleavage is right for reverse reads."""
    if not cigartuples:
        return False
    operation = cigartuples[-1][0] if is_reverse else cigartuples[0][0]
    return operation in {4, 5}


def transcript_end(read: Any) -> tuple[int, str, bool]:
    if read.is_reverse:
        return int(read.reference_end) - 1, "+", defining_end_is_clipped(read.cigartuples, True)
    return int(read.reference_start), "-", defining_end_is_clipped(read.cigartuples, False)


def eligible_mapping(read: Any, compatibility_mode: bool = False) -> tuple[bool, str]:
    if read.is_unmapped:
        return False, "unmapped"
    if read.is_secondary:
        return False, "secondary"
    if read.is_supplementary:
        return False, "supplementary"
    nh = int(read.get_tag("NH")) if read.has_tag("NH") else 0
    if not compatibility_mode and nh != 1:
        return False, "multimapping_or_missing_NH"
    return True, "eligible_primary"


def _base_rich(sequence: str, base: str, consecutive: int, window: int, minimum: int) -> bool:
    sequence = sequence.upper()
    return base * consecutive in sequence or any(
        part.count(base) >= minimum for part in (sequence[index:index + window] for index in range(max(0, len(sequence) - window + 1)))
        if len(part) == window
    )


def internal_priming_at_position(
    fasta: Any, chrom: str, position: int, strand: str,
    consecutive: int = 6, window: int = 10, minimum: int = 7,
) -> bool:
    """Test whether an observed end belongs to any qualifying A/T-rich genomic window."""
    length = int(fasta.get_reference_length(chrom))
    flank = max(consecutive, window) - 1
    start, end = max(0, position - flank), min(length, position + flank + 1)
    sequence = fasta.fetch(chrom, start, end).upper()
    local = position - start
    base = "A" if strand == "+" else "T"
    for width, required in ((consecutive, consecutive), (window, minimum)):
        first = max(0, local - width + 1)
        last = min(local, len(sequence) - width)
        for offset in range(first, last + 1):
            segment = sequence[offset:offset + width]
            if (base * consecutive in segment if width == consecutive else segment.count(base) >= required):
                return True
    return False


def centered_interval(position: int, width: int, chrom_length: int) -> tuple[int, int]:
    if width < 1:
        raise ValueError("Interval width must be positive")
    lower = (width - 1) // 2
    start = max(0, position - lower)
    end = min(chrom_length, position + (width - lower))
    return start, end


def load_rescue_sites(atlas: str | Path | None, tier: str = "core") -> dict[tuple[str, str], list[int]]:
    sites: dict[tuple[str, str], list[int]] = defaultdict(list)
    if not atlas:
        return sites
    source = Path(atlas)
    candidates: list[Path]
    if source.is_dir():
        candidates = [source / "core.bed.gz"]
        if tier == "core_plus_rescue":
            candidates.append(source / "rescue.bed.gz")
    else:
        candidates = [source]
    for path in candidates:
        if not path.is_file():
            continue
        opener = gzip.open if path.suffix == ".gz" else open
        with opener(path, "rt", encoding="utf-8") as handle:
            for line in handle:
                if not line.strip() or line.startswith(("#", "track", "browser")):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 3:
                    continue
                strand = fields[5] if len(fields) > 5 and fields[5] in {"+", "-"} else "."
                sites[(fields[0], strand)].append(int(fields[1]))
    for positions in sites.values():
        positions.sort()
    return sites


def gtf_transcript_ends(gtf: str | Path) -> dict[tuple[str, str], list[int]]:
    result: dict[tuple[str, str], list[int]] = defaultdict(list)
    path = Path(gtf)
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] not in {"transcript", "gene"} or fields[6] not in {"+", "-"}:
                continue
            start, end = int(fields[3]) - 1, int(fields[4])
            result[(fields[0], fields[6])].append(end - 1 if fields[6] == "+" else start)
    for positions in result.values():
        positions.sort()
    return result


def load_gene_models(gtf: str | Path) -> dict[str, Gene]:
    """Load genes and retain the terminal exon of every annotated transcript."""
    genes, _ = load_genes(str(gtf))
    transcript_exons: dict[tuple[str, str], list[tuple[int, int]]] = defaultdict(list)
    path = Path(gtf)
    open_gtf = gzip.open if path.suffix == ".gz" else open
    with open_gtf(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "exon":
                continue
            attributes: dict[str, str] = {}
            for item in fields[8].strip().strip(";").split(";"):
                key, _, value = item.strip().partition(" ")
                if key:
                    attributes[key] = value.strip().strip('"')
            gene_id, transcript_id = attributes.get("gene_id", ""), attributes.get("transcript_id", "")
            if gene_id in genes and transcript_id:
                transcript_exons[(gene_id, transcript_id)].append((int(fields[3]) - 1, int(fields[4])))
    terminal_by_gene: dict[str, set[tuple[int, int]]] = defaultdict(set)
    for (gene_id, _), exons in transcript_exons.items():
        gene = genes[gene_id]
        terminal = max(exons, key=lambda interval: interval[1]) if gene.strand == "+" else min(exons, key=lambda interval: interval[0])
        terminal_by_gene[gene_id].add(terminal)
    for gene_id, gene in genes.items():
        fallback = ([interval for interval in gene.exons if interval[1] == gene.end]
                    if gene.strand == "+" else [interval for interval in gene.exons if interval[0] == gene.start])
        gene.terminal_exons = sorted(terminal_by_gene.get(gene_id, set(fallback)))  # type: ignore[attr-defined]
    return genes


def rescue_overlap(
    chrom: str, strand: str, position: int, sources: Iterable[dict[tuple[str, str], list[int]]], width: int = 20,
) -> bool:
    lower = (width - 1) // 2
    upper = width - lower
    for source in sources:
        for key in ((chrom, strand), (chrom, ".")):
            positions = source.get(key, [])
            index = bisect.bisect_left(positions, position - upper + 1)
            if index < len(positions) and positions[index] - lower <= position < positions[index] + upper:
                return True
    return False


def pooled_cpm(
    sample_counts: dict[str, dict[tuple[str, str, int], int]],
) -> tuple[dict[tuple[str, str, int], float], dict[str, int]]:
    pooled: dict[tuple[str, str, int], float] = defaultdict(float)
    totals: dict[str, int] = {}
    for sample_id, counts in sample_counts.items():
        total = sum(counts.values())
        totals[sample_id] = total
        if total == 0:
            continue
        scale = 1_000_000.0 / total
        for key, count in counts.items():
            pooled[key] += count * scale
    return dict(pooled), totals


def _threshold_components(points: dict[int, float], width: int, threshold: float) -> list[tuple[int, int]]:
    events: dict[int, float] = defaultdict(float)
    for position, value in points.items():
        events[position - width + 1] += value
        events[position + 1] -= value
    retained: list[tuple[int, int]] = []
    score = 0.0
    coordinates = sorted(events)
    for index, coordinate in enumerate(coordinates[:-1]):
        score += events[coordinate]
        next_coordinate = coordinates[index + 1]
        if score > threshold and coordinate < next_coordinate:
            retained.append((coordinate, next_coordinate + width - 1))
    return merge_intervals(retained, adjacent=True)


def merge_intervals(intervals: Iterable[tuple[int, int]], adjacent: bool = True) -> list[tuple[int, int]]:
    merged: list[list[int]] = []
    for start, end in sorted(intervals):
        if end <= start:
            continue
        should_merge = bool(merged) and (start <= merged[-1][1] if adjacent else start < merged[-1][1])
        if should_merge:
            merged[-1][1] = max(merged[-1][1], end)
        else:
            merged.append([start, end])
    return [(start, end) for start, end in merged]


def _recenter(
    components: list[tuple[int, int]], points: dict[int, float], width: int, chrom_length: int,
) -> list[tuple[int, int, int]]:
    positions = sorted(points)
    result: list[tuple[int, int, int]] = []
    for start, end in components:
        left, right = bisect.bisect_left(positions, start), bisect.bisect_left(positions, end)
        members = positions[left:right]
        if not members:
            continue
        summit = min(members, key=lambda position: (-points[position], position))
        centered_start, centered_end = centered_interval(summit, width, chrom_length)
        result.append((centered_start, centered_end, summit))
    return result


def discover_sites(
    pooled: dict[tuple[str, str, int], float], chrom_lengths: dict[str, int],
    width: int = 30, threshold: float = 30.0, rounds: int = 2,
) -> list[Site]:
    if width != 30 or rounds != 2:
        raise ValueError("Mcell2019 discovery requires width=30 and rounds=2")
    grouped: dict[tuple[str, str], dict[int, float]] = defaultdict(dict)
    for (chrom, strand, position), value in pooled.items():
        grouped[(chrom, strand)][position] = grouped[(chrom, strand)].get(position, 0.0) + value
    sites: list[Site] = []
    for (chrom, strand), points in sorted(grouped.items()):
        components = _threshold_components(points, width, threshold)
        centered = _recenter(components, points, width, chrom_lengths[chrom])
        components = merge_intervals(((start, end) for start, end, _ in centered), adjacent=True)
        centered = _recenter(components, points, width, chrom_lengths[chrom])
        final_intervals = [(start, end) for start, end, _ in centered]
        if merge_intervals(final_intervals, adjacent=False) != final_intervals:
            raise RuntimeError(f"Final active-PAS intervals overlap on {chrom} {strand}")
        sites.extend(Site(chrom, strand, summit, start, end) for start, end, summit in centered)
    return sorted(sites, key=lambda site: (site.chrom, site.start, site.strand, site.summit))


def count_sites(
    sites: list[Site], sample_counts: dict[str, dict[tuple[str, str, int], int]],
) -> list[dict[str, int | str]]:
    indexed: dict[str, dict[tuple[str, str], tuple[list[int], list[int]]]] = {}
    for sample_id, counts in sample_counts.items():
        grouped: dict[tuple[str, str], list[tuple[int, int]]] = defaultdict(list)
        for (chrom, strand, position), value in counts.items():
            grouped[(chrom, strand)].append((position, value))
        sample_index: dict[tuple[str, str], tuple[list[int], list[int]]] = {}
        for key, values in grouped.items():
            values.sort()
            positions = [position for position, _ in values]
            prefix = [0]
            for _, value in values:
                prefix.append(prefix[-1] + value)
            sample_index[key] = positions, prefix
        indexed[sample_id] = sample_index
    rows: list[dict[str, int | str]] = []
    for index, site in enumerate(sites, start=1):
        row: dict[str, int | str] = {"pas_id": f"PAS_{index:08d}"}
        for sample_id, sample_index in indexed.items():
            positions, prefix = sample_index.get((site.chrom, site.strand), ([], [0]))
            left = bisect.bisect_left(positions, site.start)
            right = bisect.bisect_left(positions, site.end)
            row[sample_id] = prefix[right] - prefix[left]
        rows.append(row)
    return rows


def assign_gene(
    site: Site, genes: dict[str, Gene], extension: int = 6000,
    bins: dict[tuple[str, str, int], set[str]] | None = None, bin_size: int = 100_000,
) -> tuple[str, str, str]:
    assignments: list[tuple[str, str]] = []
    candidate_ids = (
        bins.get((site.chrom, site.strand, site.summit // bin_size), set())
        if bins is not None else genes.keys()
    )
    for gene_id in candidate_ids:
        gene = genes[gene_id]
        if gene.chrom != site.chrom or gene.strand != site.strand:
            continue
        extended_start = max(0, gene.start - extension) if gene.strand == "-" else gene.start
        extended_end = gene.end + extension if gene.strand == "+" else gene.end
        if not extended_start <= site.summit < extended_end:
            continue
        terminal = getattr(gene, "terminal_exons", None)
        if terminal is None:
            terminal = ([interval for interval in gene.exons if interval[1] == gene.end]
                        if gene.strand == "+" else [interval for interval in gene.exons if interval[0] == gene.start])
        if any(start <= site.summit < end for start, end in terminal):
            feature = "terminal_exon"
        elif any(start <= site.summit < end for start, end in gene.exons):
            feature = "other_exon"
        elif gene.start <= site.summit < gene.end:
            feature = "intron"
        else:
            feature = "downstream_extension"
        assignments.append((gene.gene_id, feature))
    if len(assignments) == 1:
        return assignments[0][0], assignments[0][1], "unique"
    if len(assignments) > 1:
        return ";".join(sorted(gene for gene, _ in assignments)), "ambiguous", "ambiguous_multi_gene"
    opposite = "-" if site.strand == "+" else "+"
    antisense_ids = bins.get((site.chrom, opposite, site.summit // bin_size), set()) if bins is not None else genes.keys()
    antisense = any(genes[gene_id].chrom == site.chrom and genes[gene_id].strand != site.strand
                    and genes[gene_id].start <= site.summit < genes[gene_id].end for gene_id in antisense_ids)
    return "", "antisense" if antisense else "intergenic", "unassigned"


def build_gene_bins(
    genes: dict[str, Gene], extension: int = 6000, bin_size: int = 100_000,
) -> dict[tuple[str, str, int], set[str]]:
    bins: dict[tuple[str, str, int], set[str]] = defaultdict(set)
    for gene in genes.values():
        start = max(0, gene.start - extension) if gene.strand == "-" else gene.start
        end = gene.end + extension if gene.strand == "+" else gene.end
        for number in range(start // bin_size, max(start, end - 1) // bin_size + 1):
            bins[(gene.chrom, gene.strand, number)].add(gene.gene_id)
    return bins


def gene_counts(
    catalog: list[dict[str, Any]], counts: list[dict[str, Any]], sample_ids: list[str],
) -> list[dict[str, Any]]:
    genes: dict[str, dict[str, int]] = defaultdict(lambda: {sample: 0 for sample in sample_ids})
    metadata = {row["pas_id"]: row for row in catalog}
    for row in counts:
        site = metadata[str(row["pas_id"])]
        if site.get("assignment_status") != "unique" or not site.get("gene_id"):
            continue
        gene_id = str(site["gene_id"])
        for sample in sample_ids:
            genes[gene_id][sample] += int(row[sample])
    return [{"gene_id": gene, **values} for gene, values in sorted(genes.items())]


def ratio_value(numerator: float, denominator: float) -> str:
    if numerator == 0 and denominator == 0:
        return "NA"
    if denominator == 0:
        return "Inf"
    return str(numerator / denominator)


def shift_direction(distal_t: float, proximal_t: float, distal_c: float, proximal_c: float) -> str:
    if (distal_t == proximal_t == 0) or (distal_c == proximal_c == 0):
        return "not_classifiable"
    left, right = distal_t * proximal_c, proximal_t * distal_c
    return "distal" if left > right else "proximal" if left < right else "no_directional_change"


def read_chrom_sizes(path: str | Path) -> dict[str, int]:
    with Path(path).open(encoding="utf-8") as handle:
        return {fields[0]: int(fields[1]) for line in handle if (fields := line.split()) and len(fields) >= 2}


def write_tsv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if fieldnames is None:
        fieldnames = list(rows[0]) if rows else []
    with path.open("w", encoding="utf-8", newline="") as handle:
        if not fieldnames:
            return
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
