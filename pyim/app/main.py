from __future__ import print_function, unicode_literals, \
                       absolute_import, division

import pandas as pd
import argparse, logging

from pyim.io.fasta import read_fasta
from pyim.io.psl import read_psl
from pyim.alignment.general import align_seqs_with_targets_par

logging.basicConfig(format='%(asctime)s \t %(message)s', level=logging.INFO)


def main():
    args = _parse_args()

    logging.info('Loading reads from fasta input file')
    fastaSeqs = list(read_fasta(args.reads))

    if args.read_limit:
        logging.info('- Limiting to %d reads' % args.read_limit)
        fastaSeqs = fastaSeqs[0:args.read_limit]

    logging.info('Loading vector definitions')
    vecSeqs = {s.seqId: s.seq for s in read_fasta(args.vectors)}

    # Align SB and T7
    logging.info('Aligning reads to SB vector')
    sbAlignments = align_seqs_with_targets_par(fastaSeqs, vecSeqs['SB'], args.num_cores, asFrame=True)

    logging.info('Aligning reads to T7 vector')
    t7Alignments = align_seqs_with_targets_par(fastaSeqs, vecSeqs['T7'], args.num_cores, asFrame=True)

    # Align barcode sequences
    logging.info('Loading barcode definitions')
    barcodeSeqs = [(s.seqId, s.seq) for s in read_fasta(args.barcodes)]

    logging.info('Aligning reads to barcodes')
    frames = []
    for bcId, bcSeq in barcodeSeqs[0:10]:
        logging.info('- Aligning reads to barcode %s' % bcId)
        bcAlignments = align_seqs_with_targets_par(fastaSeqs, bcSeq, args.num_cores, asFrame=True)
        bcAlignments['barcodeId'] = bcId
        frames.append(bcAlignments)
    bcAlignments = pd.concat(frames, axis=0)

    # Load BLAT psl alignment output
    #pslReads = read_psl(pslInput)   # Just for reference
    logging.info('Loading BLAT alignment')
    bestReads = read_psl(args.read_alignment, onlyPutative=True)

    import pdb; pdb.set_trace()


def _parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--reads', required=True)
    parser.add_argument('--vectors', required=True)
    parser.add_argument('--barcodes', required=True)
    parser.add_argument('--read_alignment', required=True)
    parser.add_argument('--num_cores', type=int, default=1)
    parser.add_argument('--read_limit', type=int, default=None)
    return parser.parse_args()


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



if __name__ == '__main__':
    main()
