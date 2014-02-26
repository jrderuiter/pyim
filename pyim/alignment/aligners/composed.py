
import logging
import pandas

from pyim.io import Sequence
from pyim.alignment.aligners.base import ReadAligner


class ChainedReadAligner(ReadAligner):

    def __init__(self, aligners, filters=None):
        super(ChainedReadAligner, self).__init__(filters)
        self.aligners = aligners
        self.logger = logging.getLogger('ChainedReadAligner')

    def _align_target(self, reads, target):
        unmapped_reads = reads

        alignments = []
        for aligner in self.aligners:
            self.logger.info('Aligning %d reads with %s' % (len(unmapped_reads), str(aligner)))
            aln_alignments, unmapped_reads = aligner.align_target(unmapped_reads, target)
            alignments.append(aln_alignments)
        self.logger.info('Leaving %d reads unaligned' % len(unmapped_reads))

        return pandas.concat(alignments, ignore_index=True), unmapped_reads


class TruncatedTargetAligner(ReadAligner):

    def __init__(self, aligner, trunc_length, filters=None):
        super(TruncatedTargetAligner, self).__init__(filters)
        self.aligner = aligner
        self.trunc_length = trunc_length

    def _align_target(self, reads, target):
        target_len = len(target)

        trunc_target = Sequence(target.name, target.seq[:self.trunc_length*-1])
        alignments, unmapped_reads = self.aligner.align_target(reads, trunc_target)

        # Fix some fields to indicate truncation has taken place
        # TODO: adjust cigar?
        type_suffix = "_truncated[%d]" % self.trunc_length
        alignments['type'] = alignments.iloc[0]['type'] + type_suffix
        alignments['identity'] = (alignments['identity'] * (target_len - self.trunc_length)) / target_len
        alignments['score'] = alignments['identity'].mul(100).round()

        return alignments, unmapped_reads

    def __str__(self):
        return '%s (%s)' % (self.__class__.__name__, self.aligner.__class__.__name__)