import unittest

from rnaends2tracks.apa_a import reverse_complement
from rnaends2tracks.mcell2019 import transcript_end


class FakeRead:
    def __init__(self, reverse, start, end, cigar):
        self.is_reverse = reverse
        self.reference_start = start
        self.reference_end = end
        self.cigartuples = cigar


class CoordinateTests(unittest.TestCase):
    def test_reverse_alignment_is_plus_transcript_rightmost_base(self):
        self.assertEqual(transcript_end(FakeRead(True, 100, 150, [(0, 50)])), (149, "+", False))

    def test_forward_alignment_is_minus_transcript_leftmost_base(self):
        self.assertEqual(transcript_end(FakeRead(False, 100, 150, [(0, 50)])), (100, "-", False))

    def test_end_defining_soft_clips(self):
        self.assertTrue(transcript_end(FakeRead(False, 100, 145, [(4, 5), (0, 45)]))[2])
        self.assertTrue(transcript_end(FakeRead(True, 100, 145, [(0, 45), (4, 5)]))[2])

    def test_hard_clips_follow_the_same_defining_side_and_other_side_is_allowed(self):
        self.assertTrue(transcript_end(FakeRead(False, 100, 145, [(5, 5), (0, 45)]))[2])
        self.assertTrue(transcript_end(FakeRead(True, 100, 145, [(0, 45), (5, 5)]))[2])
        self.assertFalse(transcript_end(FakeRead(False, 100, 145, [(0, 45), (4, 5)]))[2])
        self.assertFalse(transcript_end(FakeRead(True, 100, 145, [(4, 5), (0, 45)]))[2])

    def test_reverse_complement(self):
        self.assertEqual(reverse_complement("AACGT"), "ACGTT")


if __name__ == "__main__":
    unittest.main()
