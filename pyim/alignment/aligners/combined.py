
import pandas as pd
import numpy as np
import logging
from collections import namedtuple

from pyim.alignment.aligners.base import ReadAligner


class CombinedReadAligner(ReadAligner):

    def __init__(self, exact_aligner, inexact_aligner):
        self.exact_aligner = exact_aligner
        self.inexact_aligner = inexact_aligner

    def align_target(self, fasta_seqs, target_seq):
        return self.align_targets(fasta_seqs, [target_seq])

    def align_targets(self, fasta_seqs, target_seqs, parallel_in_targets=False, min_inexact_score=50):
        exact_frame = self._exact_alignments(fasta_seqs, target_seqs)

        unmatched_seqs = [s for s in fasta_seqs if s.seqId not in exact_frame.index]
        inexact_frame = self._inexact_alignment(unmatched_seqs, target_seqs, parallel_in_targets, select_best=True)

        return pd.concat([exact_frame, inexact_frame])

    def _exact_alignments(self, fasta_seqs, target_seqs):
        alignments = [self.exact_aligner.align_target(fasta_seqs, tseq) for tseq in target_seqs]
        return self._merge_frames(alignments)

    def _inexact_alignment(self, fasta_seqs, target_seqs, parallel_in_targets, select_best):
        logging.info('Performing inexact alignment for remaining %d reads', len(fasta_seqs))
        alignment_frame = self.inexact_aligner.align_targets(fasta_seqs, target_seqs,
                                                             parallel_in_targets=parallel_in_targets)

        if select_best:
            alignment_frame = self._select_best_inexact(alignment_frame)

        return alignment_frame

    def _select_best_inexact(self, inexact_frame):
        best_alignments = []

        grouped_inexact = inexact_frame.groupby(level=0)    # Group on index
        for name, group in grouped_inexact:
            max_score = group['score'].max()
            max_ind = np.nonzero(group['score'] == max_score)[0]

            if len(max_ind) == 1:
                best_alignments.append(group.ix[max_ind])
            else:
                logging.warning('Sequence %s has multiple alignments (%d) with same score, discarding', name, len(max_ind))

        return pd.concat(best_alignments)

    def _merge_frames(self, frames):
        frames = [frame for frame in frames if frame is not None]
        return pd.concat(frames)
