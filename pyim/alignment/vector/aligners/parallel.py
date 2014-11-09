from math import ceil
from multiprocessing import Pool
from functools import partial
from itertools import chain

from .pyim.alignment.vector.aligners.base import SequenceAligner
from pyim.alignment.general.vector.config import aligner_from_options


class ParallelSequenceAligner(SequenceAligner):

    def __init__(self, aligner, threads, options=None):
        super(ParallelSequenceAligner, self).__init__(options)

        ## If aligner is a dict it has been supplied as
        ## an option dict. Try to instantiate from options.
        if type(aligner) == dict:
            aligner = aligner_from_options(aligner)

        self.aligner = aligner
        self.threads = threads

    def __str__(self):
        return '{} ({})'.format(self.__class__.__name__, 
                                self.aligner.__class__.__name__)


class ParallelQuerySequenceAligner(ParallelSequenceAligner):

    def _align(self, queries, vector, chunk_size=None):
        ## Determine an appropriate chunk size if not specified.
        if chunk_size is None:
            chunk_size = len(queries)/float(self.threads)
            chunk_size = int(ceil(chunk_size))

        ## Perform alignment per chunk in parallel.
        pool = Pool(self.threads)
       
        func = partial(_align, vector=vector, aligner=self.aligner)
        results = pool.map(func, queries, chunk_size)
       
        pool.close()
        pool.join()

        ## Combine results into a single dictionary.
        merged = dict(chain((r.items() for r in results)))
        
        return merged


class ParallelVectorSequenceAligner(ParallelSequenceAligner):

    ## Perform alignment per vector in parallel.
    def _align_multiple(self, queries, vectors, chunk_size=1):
        pool = Pool(self.threads)
        
        func = partial(_align, queries=queries, aligner=self.aligner)
        results = pool.map(func, vectors, chunk_size)
        
        pool.close()
        pool.join()

        return results


def _align(queries, vector, aligner):
    return aligner.align(queries, vector)
