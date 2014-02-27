
import math
import functools
import multiprocessing

import pandas

from pyim.util import chunks
from pyim.alignment.vector.aligners.base import ReadAligner


class ParallelAligner(ReadAligner):

    def __init__(self, aligner, pool_size, filters=None):
        super(ParallelAligner, self).__init__(filters)
        self.aligner = aligner
        self.pool_size = pool_size

    def __str__(self):
        return '%s (%s)' % (self.__class__.__name__, self.aligner.__class__.__name__)


class ParallelReadAligner(ParallelAligner):

    def _align_target(self, reads, target):
        chunk_size = int(math.ceil(len(reads)/float(self.pool_size)))
        read_chunks = chunks(reads, chunk_size)

        # Do alignment and retrieve results
        pool = multiprocessing.Pool(self.pool_size)
        func = functools.partial(_align_target, target_seq=target, aligner=self.aligner)
        results = pool.map(func, read_chunks, 1)
        pool.close()
        pool.join()

        # Concatenate results, filter nones
        alignments, unmapped = zip(*results)

        # First the alignments
        alns_not_none = [aln for aln in alignments if aln is not None]
        if alns_not_none:
            alns_concat = pandas.concat(alns_not_none, ignore_index=True)
        else:
            alns_concat = None

        # Then the unmapped reads
        unmapped_concat = [item for sublist in unmapped for item in sublist]

        return alns_concat, unmapped_concat


class ParallelTargetAligner(ParallelAligner):

    def align_targets(self, reads, targets):
        if type(targets) == dict:
            targets = targets.values()

        chunk_size = int(math.ceil(len(targets)/float(self.pool_size)))

        pool = multiprocessing.Pool(self.pool_size)
        func = functools.partial(_align_target_tup, reads=reads, aligner=self.aligner)
        tup_results = pool.map(func, targets, chunk_size)
        pool.close()
        pool.join()

        names, results = zip(*tup_results)
        alignments, unmapped = zip(*results)

        alignment_dict = dict(zip(names, alignments))
        unmapped_dict = dict(zip(names, unmapped))

        return alignment_dict, unmapped_dict


def _align_target(reads, target_seq, aligner):
    return aligner.align_target(reads, target_seq)


def _align_target_tup(target_seq, reads, aligner):
    result = _align_target(reads, target_seq, aligner)
    return target_seq.name, result

