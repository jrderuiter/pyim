
import multiprocessing, functools

from pyim.alignment.base import Alignment
from pyim.alignment.algorithms.waterman import water

def align_seqs_with_targets(fastaSeqs, targetSeq, asFrame=False):
    alignments = {fSeq.seqId: water(fSeq.seq, targetSeq) for fSeq in fastaSeqs}
    if asFrame:
        alignments = Alignment.as_dataframe(alignments)
    return alignments


def align_seqs_with_targets_par(fastaSeqs, targetSeq, numProcesses, asFrame=False):
    water_partial = functools.partial(water_tuple, targetSeq=targetSeq)

    pool = multiprocessing.Pool(numProcesses)
    rawAlignments = pool.map(water_partial, fastaSeqs)

    alignments = dict(rawAlignments)
    if asFrame:
        alignments = Alignment.as_dataframe(alignments)

    return alignments


def water_tuple(fSeq, targetSeq):
    return fSeq.seqId, water(fSeq.seq, targetSeq)
