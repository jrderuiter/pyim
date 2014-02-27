import unittest

from pyim.alignment.vector.algorithms.waterman import water


class TestWaterman(unittest.TestCase):

    EXO_MATCH_SCORE = 5
    EXO_MISMATCH_SCORE = -4
    EXO_GAP_SCORE = -12


    def test_ungapped(self):
        readSeq = 'ACGTCTCATCGTGTATGTAACTTCCGACTTCAACTGTATAGGGAGTCCTACTAGTGAGTC'
        sbSeq = 'GTGTATGTAAACTTCCGACTTCAAC'

        ja, jb, n = water(readSeq, sbSeq, gapScore=-100)

        self.assertEqual(ja, 18)
        self.assertEqual(jb, 9)
        self.assertEqual(n, 16)
        self.assertEqual(readSeq[ja:ja + n], sbSeq[jb:jb + n])
        self.assertEqual(readSeq[ja:ja + n], 'AACTTCCGACTTCAAC')

    def test_nohit(self):
        readSeq = 'TTTTTTTTT'
        sbSeq = 'AA'

        ja, jb, n = water(readSeq, sbSeq)

        self.assertEqual(ja, 0)
        self.assertEqual(jb, 0)
        self.assertEqual(n, 0)

    def test_mismatch_exonerate_one(self):
        readSeq = 'ACGCTCGACAGTGTATGAAAACTTCCGACTTCAACTGTATAGGGATCCTCTAGCCTATAG'
        sb010Seq = 'ACGCTCGACAGTGTATGTAA'

        ja, jb, n = water(readSeq, sb010Seq, matchScore=self.EXO_MATCH_SCORE,
                          mismatchScore=self.EXO_MISMATCH_SCORE, gapScore=self.EXO_GAP_SCORE)
        self.assertEqual(readSeq[ja:ja + n], 'ACGCTCGACAGTGTATGAAA')
        self.assertEqual(sb010Seq[jb:jb + n], 'ACGCTCGACAGTGTATGTAA')

    def test_mismatch_exonerate_two(self):
        readSeq = 'ACTGTACAGTGTGTATGTAAACTTCCGACTTCAACTGTATCCTATAGTGAGTCGTATTA'
        sb010Seq = 'CAGGTCCAGTGTGTATGTAA'

        ja, jb, n = water(readSeq, sb010Seq, matchScore=self.EXO_MATCH_SCORE,
                          mismatchScore=self.EXO_MISMATCH_SCORE, gapScore=self.EXO_GAP_SCORE)
        self.assertEqual(readSeq[ja:ja + n], 'GTACAGTGTGTATGTAA')
        self.assertEqual(sb010Seq[jb:jb + n], 'GTCCAGTGTGTATGTAA')


