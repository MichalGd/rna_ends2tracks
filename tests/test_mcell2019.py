import unittest

from rnaends2tracks.apa_a import Gene
from rnaends2tracks.mcell2019 import (
    Site,
    assign_gene,
    build_gene_bins,
    centered_interval,
    count_sites,
    discover_sites,
    eligible_mapping,
    gene_counts,
    internal_priming_at_position,
    pooled_cpm,
    ratio_value,
    rescue_overlap,
    shift_direction,
)


class FakeFasta:
    def __init__(self, sequence): self.sequence = sequence
    def get_reference_length(self, _chrom): return len(self.sequence)
    def fetch(self, _chrom, start, end): return self.sequence[start:end]


class FakeRead:
    def __init__(self, *, unmapped=False, secondary=False, supplementary=False, duplicate=False, nh=1):
        self.is_unmapped = unmapped
        self.is_secondary = secondary
        self.is_supplementary = supplementary
        self.is_duplicate = duplicate
        self._nh = nh

    def has_tag(self, tag): return tag == "NH" and self._nh is not None
    def get_tag(self, _tag): return self._nh


class Mcell2019Tests(unittest.TestCase):
    def test_even_width_interval_is_documented_bed_convention(self):
        self.assertEqual(centered_interval(100, 30, 1000), (86, 116))
        self.assertEqual(centered_interval(2, 30, 1000), (0, 18))
        self.assertEqual(centered_interval(998, 30, 1000), (984, 1000))

    def test_strict_threshold_and_low_coordinate_tie_break(self):
        equal = {("chr1", "+", 100): 31.0, ("chr1", "+", 101): 31.0}
        sites = discover_sites(equal, {"chr1": 1000})
        self.assertEqual((sites[0].summit, sites[0].start, sites[0].end), (100, 86, 116))
        self.assertEqual(discover_sites({("chr1", "+", 100): 30.0}, {"chr1": 1000}), [])
        minus = discover_sites({("chr1", "-", 100): 31.0, ("chr1", "-", 101): 31.0}, {"chr1": 1000})
        self.assertEqual(minus[0].summit, 100)

    def test_cpm_is_sample_level_across_both_strands(self):
        pooled, totals = pooled_cpm({"S1": {("chr1", "+", 1): 1, ("chr1", "-", 2): 1}})
        self.assertEqual(totals["S1"], 2)
        self.assertEqual(pooled[("chr1", "+", 1)], 500000)
        self.assertEqual(pooled[("chr1", "-", 2)], 500000)

    def test_internal_priming_windows_include_observed_position(self):
        fasta = FakeFasta("CCCCCAAAAAACCCCC")
        self.assertTrue(internal_priming_at_position(fasta, "chr1", 7, "+"))
        self.assertFalse(internal_priming_at_position(fasta, "chr1", 1, "+"))
        seven_of_ten = FakeFasta("CCCAAACAACAAACCCC")
        self.assertTrue(internal_priming_at_position(seven_of_ten, "chr1", 8, "+"))
        minus = FakeFasta("CCCTTTCTTCTTTCCCC")
        self.assertTrue(internal_priming_at_position(minus, "chr1", 8, "-"))

    def test_mapping_policy_excludes_nonunique_and_retains_duplicates(self):
        self.assertEqual(eligible_mapping(FakeRead(duplicate=True)), (True, "eligible_primary"))
        self.assertEqual(eligible_mapping(FakeRead(nh=2)), (False, "multimapping_or_missing_NH"))
        self.assertEqual(eligible_mapping(FakeRead(secondary=True)), (False, "secondary"))
        self.assertEqual(eligible_mapping(FakeRead(supplementary=True)), (False, "supplementary"))
        self.assertEqual(eligible_mapping(FakeRead(unmapped=True)), (False, "unmapped"))

    def test_rescue_interval_and_gene_assignment_contracts(self):
        source = {("chr1", "+"): [100]}
        self.assertTrue(rescue_overlap("chr1", "+", 91, (source,), 20))
        self.assertTrue(rescue_overlap("chr1", "+", 110, (source,), 20))
        self.assertFalse(rescue_overlap("chr1", "+", 111, (source,), 20))
        plus = Gene("plus", "chr1", 100, 200, "+", exons=[(100, 130), (170, 200)])
        plus.terminal_exons = [(170, 200)]
        minus = Gene("minus", "chr1", 300, 400, "-", exons=[(300, 330), (370, 400)])
        minus.terminal_exons = [(300, 330)]
        genes = {gene.gene_id: gene for gene in (plus, minus)}
        bins = build_gene_bins(genes, extension=6000)
        self.assertEqual(assign_gene(Site("chr1", "+", 250, 236, 266), genes, 6000, bins),
                         ("plus", "downstream_extension", "unique"))
        self.assertEqual(assign_gene(Site("chr1", "-", 250, 236, 266), genes, 6000, bins),
                         ("minus", "downstream_extension", "unique"))
        self.assertEqual(assign_gene(Site("chr1", "+", 150, 136, 166), genes, 6000, bins),
                         ("plus", "intron", "unique"))

    def test_ambiguous_sites_are_excluded_from_c4(self):
        catalog = [
            {"pas_id": "P1", "gene_id": "G1", "assignment_status": "unique"},
            {"pas_id": "P2", "gene_id": "G1;G2", "assignment_status": "ambiguous_multi_gene"},
        ]
        counts = [{"pas_id": "P1", "S1": 5}, {"pas_id": "P2", "S1": 9}]
        self.assertEqual(gene_counts(catalog, counts, ["S1"]), [{"gene_id": "G1", "S1": 5}])

    def test_each_end_is_counted_in_at_most_one_nonoverlapping_site(self):
        sites = discover_sites({("chr1", "+", 100): 40, ("chr1", "+", 200): 40}, {"chr1": 1000})
        rows = count_sites(sites, {"S1": {("chr1", "+", 100): 3, ("chr1", "+", 200): 4}})
        self.assertEqual(sum(int(row["S1"]) for row in rows), 7)

    def test_zero_safe_ratios_and_direction(self):
        self.assertEqual(ratio_value(0, 0), "NA")
        self.assertEqual(ratio_value(1, 0), "Inf")
        self.assertEqual(shift_direction(2, 1, 1, 2), "distal")


if __name__ == "__main__":
    unittest.main()
