
from __future__ import print_function, absolute_import, division

import multiprocessing
import functools
import sys
import pandas as pd

from pyim.alignment.aligners.base import ReadAligner

PARALLEL_CHUNK_SIZE = 200


class ParallelReadAligner(ReadAligner):

    def __init__(self, num_processes=1, align_func=None):
        if align_func is None:
            align_func = _not_implemented

        self.align_func = align_func
        self.num_processes = num_processes

    def align_target(self, fasta_seqs, target_seq):
        if self.num_processes == 1:
            alignments = _sequential(fasta_seqs, target_seq, self.align_func)
        else:
            alignments = self._parallel_in_reads_progress(fasta_seqs, target_seq)
        return self._to_frame(alignments)

    def align_targets(self, fasta_seqs, target_seqs, parallel_in_targets=False):
        if parallel_in_targets:
            alignments = self._parallel_in_targets(fasta_seqs, target_seqs)
        else:
            alignments = super(ParallelReadAligner, self).align_targets(fasta_seqs, target_seqs)
        return alignments

    def _parallel_in_reads(self, fasta_seqs, target_seq):
        aln_partial = functools.partial(self.align_func, target=target_seq)

        pool = multiprocessing.Pool(self.num_processes)
        alignments = pool.map(aln_partial, fasta_seqs)
        pool.close()
        pool.join()

        return alignments

    def _parallel_in_reads_progress(self, fasta_seqs, target_seq):
        num_reads = float(len(fasta_seqs))
        aln_partial = functools.partial(self.align_func, target=target_seq)

        pool = multiprocessing.Pool(self.num_processes)

        alignments = []
        _print_unbuffered('\r\tProcessed 0%% of %d reads' % num_reads, end='')
        for i, aln in enumerate(pool.imap_unordered(aln_partial, fasta_seqs, PARALLEL_CHUNK_SIZE), 1):
            alignments.append(aln)
            if i % (PARALLEL_CHUNK_SIZE * self.num_processes) == 0:
                _print_unbuffered('\r\tProcessed %3.2f%% of %d reads' % ((i / num_reads) * 100, num_reads), end='')
        _print_unbuffered('\r\tProcessed 100%% of %d reads' % num_reads)

        pool.close()
        pool.join()

        return alignments

    def _parallel_in_targets(self, fasta_seqs, target_seqs):
        seq_partial = functools.partial(_sequential_targets, fasta_seqs=fasta_seqs, align_func=self.align_func)
        num_targets = len(target_seqs)

        pool = multiprocessing.Pool(self.num_processes)

        alignments = []
        _print_unbuffered('\r\tProcessed 0%% of %d target sequences' % num_targets, end='')
        for i, aln in enumerate(pool.imap_unordered(seq_partial, target_seqs, 1), 1):
            alignments.append(aln)
            if i % (PARALLEL_CHUNK_SIZE * self.num_processes) == 0:
                _print_unbuffered('\r\tProcessed %3.2f%% of %d target sequences' % ((i / num_targets) * 100, num_targets), end='')
        _print_unbuffered('\r\tProcessed 100%% of %d target sequences' % num_targets)

        pool.close()
        pool.join()

        return pd.concat([self._to_frame(aln) for aln in alignments])


def _sequential(fasta_seqs, target_seq, align_func):
    return [align_func(fseq, target_seq) for fseq in fasta_seqs]

def _sequential_targets(target_seq, fasta_seqs, align_func):
    return [align_func(fseq, target_seq) for fseq in fasta_seqs]

def _print_unbuffered(message, end='\n'):
    print(message, end=end)
    sys.stdout.flush()

def _not_implemented(readSeq, targetSeq):
    raise NotImplementedError
