
import unittest
from pyim.alignment.algorithms.lcs import longest_common_subsequence


class TestLCS(unittest.TestCase):

    def test_submatch(self):
        readSeq = 'ACGTCTCATCGTGTATGTAACTTCCGACTTCAACTGTATAGGGAGTCCTACTAGTGAGTC'
        sbSeq = 'GTGTATGTAAACTTCCGACTTCAAC'

        ja, jb, n = longest_common_subsequence(readSeq, sbSeq)

        self.assertEqual(ja, 18)
        self.assertEqual(jb, 9)
        self.assertEqual(n, 16)
        self.assertEqual(readSeq[ja:ja + n], sbSeq[jb:jb + n])
        self.assertEqual(readSeq[ja:ja + n], 'AACTTCCGACTTCAAC')

    def test_full_match(self):
        readSeq = 'ACGCTCGACAGTGTATGTAACTTCGACTTCAACGCCTATAGTGAGTCGTATTA'
        t7Seq = 'CCTATAGTGAGTCGTATTA'

        ja, jb, n = longest_common_subsequence(readSeq, t7Seq)

        self.assertEqual(ja, 34)
        self.assertEqual(jb, 0)
        self.assertEqual(n, len(t7Seq))
        self.assertEqual(readSeq[ja:ja + n], 'CCTATAGTGAGTCGTATTA')
        self.assertEqual(readSeq[ja:ja + n], t7Seq[jb:jb + n])

    def test_custom_match(self):
        readSeq = 'ACGCCCTATAGTGAGTCGTANTATATTA'
        t7Seq = 'CCTATAGTGAGTCGTATTA'

        ja, jb, n = longest_common_subsequence(readSeq, t7Seq)
        self.assertEqual(readSeq[ja:ja + n], 'CCTATAGTGAGTCGTA')

        def match(x,y):
            if x == 'N' or y == 'N':
                return True
            return x == y

        ja, jb, n = longest_common_subsequence(readSeq, t7Seq, matchFunc=match)
        self.assertEqual(readSeq[ja:ja + n], 'CCTATAGTGAGTCGTANTA')

