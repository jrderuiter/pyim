from __future__ import print_function, unicode_literals, \
                       absolute_import, division
import os.path as path
import multiprocessing
import functools
import logging

import pandas as pd

from pyim.io.fasta import read_fasta
from pyim.io.psl import read_psl
from pyim.alignment.base import Alignment
from pyim.alignment.algorithms.waterman import water
logging.basicConfig(format='%(asctime)s \t %(message)s', level=logging.INFO)

DATA_DIR = '/Users/j.d.ruiter/Dropbox/Projects/NKI/im/PyIM-pipeline/data'
BASE_DIR = '/Volumes/Datastore/Julian/Datasets/sjors-im/Pool3'

fastaInput = path.join(BASE_DIR, 'Pool3.1.TCA.454Reads.fna')
pslInput = path.join(BASE_DIR, 'POOL3.1.psl')
vectorDefFile = path.join(DATA_DIR, 'vec.fa')


def main():
    logging.info('Loading reads from fasta input file')
    fastaSeqs = list(read_fasta(fastaInput))
    logging.info('Loading vector definitions')
    vecSeqs = {s.seqId: s.seq for s in read_fasta(vectorDefFile)}

    # Align SB and T7
    logging.info('Aligning reads to SB vector')
    sbAlignments = align_seqs_with_targets_par(fastaSeqs, vecSeqs['SB'], asFrame=True)

    import pdb; pdb.set_trace()
    logging.info('Aligning reads to T7 vector')
    t7Alignments = align_seqs_with_targets(fastaSeqs, vecSeqs['T7'], asFrame=True)

    # Align barcode sequences
    logging.info('Loading barcode definitions')
    barcodeSeqs = [(s.seqId, s.seq) for s in read_fasta(path.join(DATA_DIR, 'SBbarcodes.fa'))]

    logging.info('Aligning reads to barcodes')
    frames = []
    for bcId, bcSeq in barcodeSeqs[0:10]:
        logging.info('- Aligning reads to barcode %s' % bcId)
        bcAlignments = align_seqs_with_targets(fastaSeqs, bcSeq, asFrame=True)
        bcAlignments['barcodeId'] = bcId
        frames.append(bcAlignments)
    bcAlignments = pd.concat(frames, axis=0)

    # Load BLAT psl alignment output
    #pslReads = read_psl(pslInput)   # Just for reference
    logging.info('Loading BLAT alignment')
    bestReads = read_psl(pslInput, onlyPutative=True)

    pdb.set_trace()


def print_vector_alignments(readSeq, vecAlignments):
    print(readSeq)
    readSeq = list(readSeq)
    for vecName, vecAln in vecAlignments.items():
        start, end = vecAln['query_start'], vecAln['query_end']

        size = end - start
        paddingSize = size - len(vecName) - 3
        vecStr = '[-' + vecName + ('-' * paddingSize) + ']'
        readSeq[start:end] = list(vecStr)
    print(''.join(readSeq))


def align_seqs_with_targets(fastaSeqs, targetSeq, asFrame=False):
    alignments = {fSeq.seqId: water(fSeq.seq, targetSeq) for fSeq in fastaSeqs}
    if asFrame:
        alignments = Alignment.as_dataframe(alignments)
    return alignments


def align_seqs_with_targets_par(fastaSeqs, targetSeq, asFrame=False):
    water_partial = functools.partial(water_tuple, targetSeq=targetSeq)

    pool = multiprocessing.Pool(4)
    rawAlignments = pool.map(water_partial, fastaSeqs)

    alignments = dict(rawAlignments)
    if asFrame:
        alignments = Alignment.as_dataframe(alignments)

    return alignments


def water_tuple(fSeq, targetSeq):
    return fSeq.seqId, water(fSeq.seq, targetSeq)


if __name__ == '__main__':
    main()
