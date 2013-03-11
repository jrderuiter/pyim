#!/usr/bin/env python

from __future__ import print_function, unicode_literals, \
                       absolute_import, division

import pandas as pd
import argparse, logging

from pyim.io.fasta import read_fasta
from pyim.io.psl import read_psl
from pyim.alignment.aligner import WatermanAligner
from pyim.alignment.cache import AlignmentCache

logging.basicConfig(format='[%(asctime)s] %(levelname)-8s %(message)s', level=logging.INFO, datefmt='%H:%M:%S')


def main(fasta_file, vector_file, barcode_file, num_cores, read_limit):
    cache = AlignmentCache(fasta_file)
    aligner = WatermanAligner(num_processes=num_cores, as_frame=True)

    fastaSeqs = _load_fasta(fasta_file, read_limit)

    logging.info('Loading vector definitions')
    vecSeqs = {s.seqId: s.seq for s in read_fasta(vector_file)}

    # Align SB and T7
    sbAlignments = _align_reads(fastaSeqs, 'SB', vecSeqs['SB'], aligner, cache)
    t7Alignments = _align_reads(fastaSeqs, 'T7', vecSeqs['T7'], aligner, cache)

    # Align barcode sequences
    logging.info('Loading barcode definitions')
    barcodeSeqs = [(s.seqId, s.seq) for s in read_fasta(barcode_file)]

    logging.info('Aligning reads to barcodes')
    frames = []
    for bcId, bcSeq in barcodeSeqs[0:15]:
        bcAlignments = _align_reads(fastaSeqs, bcId, bcSeq, aligner, cache)
        bcAlignments['barcodeId'] = bcId
        frames.append(bcAlignments)
    #bcAlignments = pd.concat(frames, axis=0)

    # Load BLAT psl alignment output
    #pslReads = read_psl(pslInput)   # Just for reference
    #logging.info('Loading BLAT alignment')
    #bestReads = read_psl(args.read_alignment, onlyPutative=True)


def _parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--reads', required=True)
    parser.add_argument('--vectors', required=True)
    parser.add_argument('--barcodes', required=True)
    parser.add_argument('--read_alignment', required=True)
    parser.add_argument('--num_cores', type=int, default=1)
    parser.add_argument('--read_limit', type=int, default=None)
    return parser.parse_args()


def _load_fasta(fasta_file, read_limit):
    logging.info('Loading reads from fasta input file')
    fastaSeqs = list(read_fasta(fasta_file))

    if read_limit:
        logging.info('- Limiting to %d reads' % read_limit)
        fastaSeqs = fastaSeqs[0:read_limit]

    return fastaSeqs


def _filter_input_reads(fasta_seqs, sb_alignments):
    starts_at_begin = sb_alignments['target_start'] == 0
    match_long_enough = sb_alignments['target_end'] >= 23
    alignment_sufficient = sb_alignments['score'] >= 75
    good_match = sb_alignments['identity'] >= 95

    seq_dict = {s.seqId: s for s in fasta_seqs}
    valid_mask = starts_at_begin & match_long_enough & \
                 alignment_sufficient & good_match
    valid_ids = list(valid_mask.index[valid_mask])
    valid_reads = [seq_dict[read_id] for read_id in valid_ids]

    return valid_reads


def _align_reads(fastaSeqs, vecName, vecSeq, aligner, cache):
    if vecName in cache:
        logging.info('Loaded alignments from cache for %s vector' % vecName)
        alignments = cache.get(vecName)
    else:
        logging.info('Aligning reads to %s vector' % vecName)
        alignments = aligner.align_target(fastaSeqs, vecSeq)
        cache.add(vecName, alignments)
    return alignments


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
    args = _parse_args()
    main(args.reads, args.vectors, args.barcodes, args.num_cores, args.read_limit)
