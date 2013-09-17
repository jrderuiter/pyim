

import pandas

from pyim.alignment.aligners.model import Alignment, alignments_to_frame
from pyim.alignment.aligners.algorithms.waterman import water
from pyim.alignment.aligners.algorithms.lcs import longest_common_subsequence


class ReadAligner(object):

    def align_target(self, fasta_seqs, target_seq):
        raise NotImplementedError

    def align_targets(self, fasta_seqs, target_seqs):
        alignments, unmapped_reads = {}, {}

        for target_seq in target_seqs:
            tgt_alignments, tgt_unmapped = self.align_target(fasta_seqs, target_seq)
            alignments[target_seq.name] = tgt_alignments
            unmapped_reads[target_seq.name] = tgt_unmapped

        return alignments, unmapped_reads


class ExactReadAligner(ReadAligner):

    def align_target(self, fasta_seqs, target):
        tgt_seq, tgt_len = target.seq, len(target)

        alignments, unmapped_seqs = {}, []
        for read in fasta_seqs:
            if tgt_seq in read.seq:
                start = read.seq.index(tgt_seq)
                cigar_str = '%dM' % tgt_len

                alignment = Alignment(read.name, start, start + tgt_len, read.seq,
                                      target.name, 0, tgt_len, tgt_seq,
                                      100, 1.0, cigar_str, 'exact')

                alignments[read.name] = alignment
            else:
                unmapped_seqs.append(read)

        alignment_frame = alignments_to_frame(alignments)
        return alignment_frame, unmapped_seqs


class InexactReadAligner(ReadAligner):

    def __init__(self, min_score=0):
        super(InexactReadAligner, self).__init__()
        self.min_score = min_score

    def _align(self, fasta_seq, target_seq):
        raise NotImplementedError

    def align_target(self, fasta_seqs, target_seq):
        alignments, unmapped_reads = {}, []

        for fasta_seq in fasta_seqs:
            alignment = self._align(fasta_seq, target_seq)

            if alignment.score >= self.min_score:
                alignments[fasta_seq.name] = alignment
            else:
                unmapped_reads.append(fasta_seq)

        alignment_frame = alignments_to_frame(alignments)
        return alignment_frame, unmapped_reads


class WatermanAligner(InexactReadAligner):

    def __init__(self, min_score=0, match_score=5, mismatch_score=-9, gap_score=-5, verbose=False):
        super(WatermanAligner, self).__init__(min_score)
        self.verbose = verbose

        self.match_score = match_score
        self.mismatch_score = mismatch_score
        self.gap_score = gap_score

    def _align(self, fasta_seq, target_seq):
        return water(fasta_seq, target_seq, matchScore=self.match_score,
                     mismatchScore=self.mismatch_score, gapScore=self.gap_score,
                     verbose=self.verbose)


class LCSAligner(InexactReadAligner):

    def _align(self, fasta_seq, target_seq):
        return longest_common_subsequence(fasta_seq, target_seq)



