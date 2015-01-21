from math import ceil
from functools import partial
import multiprocessing

from pyim.common.alignment.vector.aligners import SequenceAligner
from pyim.common.util import chunks


class ParallelSequenceAligner(SequenceAligner):

    def __init__(self, aligner, threads):
        super(ParallelSequenceAligner, self).__init__()

        if isinstance(aligner, ParallelSequenceAligner):
            raise ValueError(('Cannot parallelize an aligner that is already '
                              'operating in parallel ({})').
                             format(aligner.__class__.__name__))

        self.aligner = aligner
        self.threads = threads

    def __str__(self):
        return '{} ({})'.format(self.__class__.__name__, 
                                self.aligner.__class__.__name__)


class ParallelQuerySequenceAligner(ParallelSequenceAligner):

    def _align(self, queries, vector, chunk_size=None):
        # Determine an appropriate chunk size if not specified.
        if chunk_size is None:
            chunk_size = len(queries)/float(self.threads)
            chunk_size = int(ceil(chunk_size))

        # Divide queries into chunks.
        query_chunks = chunks(queries, chunk_size)

        # Perform alignment per chunk in parallel.
        pool = multiprocessing.Pool(self.threads)
       
        func = partial(_align, vector=vector, aligner=self.aligner)
        results = pool.map(func, query_chunks)
       
        pool.close()
        pool.join()

        # Combine results into a single dictionary.
        merged = {k: v for d in results for k, v in d.items()}

        return merged


class ParallelVectorSequenceAligner(ParallelSequenceAligner):

    # Perform alignment per vector in parallel.
    def _align_multiple(self, queries, vectors):
        pool = multiprocessing.Pool(self.threads)
        
        func = partial(_align, queries=queries, aligner=self.aligner)
        results = pool.map(func, vectors, 1)
        
        pool.close()
        pool.join()

        return results


def _align(queries, vector, aligner):
    return aligner.align(queries, vector)
