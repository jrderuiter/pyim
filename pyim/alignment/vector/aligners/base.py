import pandas

from pyim.alignment.vector.model import Alignment
from pyim.alignment.vector.algorithms.waterman import water
from pyim.alignment.vector.algorithms.lcs import longest_common_subsequence


class ReadAligner(object):

    def __init__(self, filters=None):
        self.filters = filters

    def align_target(self, reads, target):
        alignment, unmapped = self._align_target(reads, target)

        if type(alignment) == list or type(alignment) == dict:
            alignment = self._alignments_to_frame(alignment)

        if self.filters is not None:
            alignment, unmapped = self._apply_filters(alignment, unmapped)

        return alignment, unmapped

    def _align_target(self, reads, target):
        raise NotImplementedError

    def align_targets(self, reads, targets):
        if type(targets) == dict:
            targets = targets.values()

        alignments, unmapped_reads = {}, {}
        for target in targets:
            tgt_alignments, tgt_unmapped = self.align_target(reads, target)
            alignments[target.name] = tgt_alignments
            unmapped_reads[target.name] = tgt_unmapped

        return alignments, unmapped_reads

    def _alignments_to_frame(self, alns):
        if len(alns) == 0:
            return None
        values = list(alns.values()) if type(alns) == dict else alns
        return pandas.DataFrame(values, columns=values[0]._fields)

    def _apply_filters(self, alignment, unmapped):
        filtered_alignment = alignment
        for filt in self.filters:
            filtered_alignment, _, dropped = filt.apply(filtered_alignment)
            unmapped = unmapped + dropped
        return filtered_alignment, unmapped

    def __str__(self):
        return self.__class__.__name__


class ExactReadAligner(ReadAligner):

    def _align_target(self, reads, target):
        tgt_seq, tgt_len = target.seq, len(target.seq)

        alignments, unmapped_seqs = {}, []
        for read in reads:
            if tgt_seq in read.seq:
                start = read.seq.index(tgt_seq)
                cigar_str = '%dM' % tgt_len

                alignment = Alignment(read.name, start, start + tgt_len, read.seq,
                                      target.name, 0, tgt_len, tgt_seq,
                                      100, 1.0, cigar_str, 'exact')
                alignments[read.name] = alignment
            else:
                unmapped_seqs.append(read)

        return alignments, unmapped_seqs


class InexactReadAligner(ReadAligner):

    def __init__(self, filters=None, min_score=0):
        super(InexactReadAligner, self).__init__(filters)
        self.min_score = min_score

    def _align(self, fasta_seq, target_seq):
        raise NotImplementedError

    def _align_target(self, fasta_seqs, target_seq):
        alignments, unmapped_reads = {}, []

        for fasta_seq in fasta_seqs:
            alignment = self._align(fasta_seq, target_seq)

            if alignment.score >= self.min_score:
                alignments[fasta_seq.name] = alignment
            else:
                unmapped_reads.append(fasta_seq)

        return alignments, unmapped_reads


class WatermanAligner(InexactReadAligner):

    def __init__(self, filters=None, min_score=0, match_score=5, mismatch_score=-9, gap_score=-5, verbose=False):
        super(WatermanAligner, self).__init__(filters, min_score)

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



