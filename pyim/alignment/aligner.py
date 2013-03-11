
from __future__ import print_function, unicode_literals, \
    absolute_import, division

import multiprocessing, functools

from pyim.alignment.base import Alignment
from pyim.alignment.algorithms.waterman import water
from pyim.alignment.algorithms.lcs import longest_common_subsequence


class ReadAligner(object):

    def __init__(self, num_processes=1, as_frame=False, align_func=None):
        if align_func is None:
            align_func = _not_implemented

        self.align_func = align_func
        self.num_processes = num_processes
        self.as_frame = as_frame

    def align_target(self, fastaSeqs, targetSeq, as_frame=None):
        if self.num_processes == 1:
            alignments = self._sequential(fastaSeqs, targetSeq)
        else:
            alignments = self._parallel_in_reads_progress(fastaSeqs, targetSeq, self.num_processes)

        as_frame = as_frame if as_frame is not None else self.as_frame
        if as_frame:
            alignments = Alignment.as_dataframe(alignments)

        return alignments

    def _align_tuple(self, readSeq, targetSeq):
        return (readSeq.seqId, self.align_func(readSeq.seq, targetSeq))

    def _sequential(self, fastaSeqs, targetSeq):
        return {fSeq.seqId: self.align_func(fSeq.seq, targetSeq) for fSeq in fastaSeqs}

    def _parallel_in_reads(self, fastaSeqs, targetSeq, numProcesses):
        aln_partial = functools.partial(_align_tuple, targetSeq=targetSeq, alignFunc=self.align_func)

        pool = multiprocessing.Pool(numProcesses)
        rawAlignments = pool.map(aln_partial, fastaSeqs)

        return dict(rawAlignments)

    def _parallel_in_reads_progress(self, fastaSeqs, targetSeq, numProcesses):
        num_reads = float(len(fastaSeqs))
        aln_partial = functools.partial(_align_tuple, targetSeq=targetSeq, alignFunc=self.align_func)

        pool = multiprocessing.Pool(numProcesses)

        rawAlignments = []
        for i, aln in enumerate(pool.imap_unordered(aln_partial, fastaSeqs), 1):
            rawAlignments.append(aln)
            if i % 500 == 0:
                print('\r\tProcessed %3.2f%% of %d reads' % ((i / num_reads) * 100, num_reads), end='')
        print('')

        return dict(rawAlignments)


def _not_implemented(readSeq, targetSeq):
    raise NotImplementedError

def _align_tuple(readSeq, targetSeq, alignFunc):
    return (readSeq.seqId, alignFunc(readSeq.seq, targetSeq))


class WatermanAligner(ReadAligner):

    def __init__(self, num_processes=1, as_frame=False):
        super(WatermanAligner, self).__init__(num_processes, as_frame)
        self.align_func = water


class LCSAligner(ReadAligner):

    def __init__(self, num_processes=1, as_frame=False):
        super(LCSAligner, self).__init__(num_processes, as_frame)
        self.align_func = longest_common_subsequence