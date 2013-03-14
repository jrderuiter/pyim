
from pyim.alignment.aligners.parallel import ParallelReadAligner
from pyim.alignment.algorithms.waterman import water
from pyim.alignment.algorithms.lcs import longest_common_subsequence


class InexactParallelReadAligner(ParallelReadAligner):

    def __init__(self, num_processes=1, align_func=None, min_score=None):
        super(InexactParallelReadAligner, self).__init__(num_processes=num_processes, align_func=align_func)
        self.min_score = min_score

    def align_target(self, fasta_seqs, target_seq):
        alignments = super(InexactParallelReadAligner, self).align_target(fasta_seqs, target_seq)
        if self.min_score is not None:
            alignments = alignments[alignments['score'] >= self.min_score]
        return alignments


class WatermanAligner(InexactParallelReadAligner):

    def __init__(self, num_processes=1, min_score=None):
        super(WatermanAligner, self).__init__(num_processes=num_processes, align_func=water,
                                              min_score=min_score)


class LCSAligner(InexactParallelReadAligner):

    def __init__(self, num_processes=1, min_score=None):
        super(LCSAligner, self).__init__(num_processes=num_processes, align_func=longest_common_subsequence,
                                         min_score=min_score)