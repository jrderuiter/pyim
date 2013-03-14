#!/usr/bin/env python

from __future__ import print_function, unicode_literals, \
                       absolute_import, division
import argparse
import logging
import pandas as pd
import numpy as np

from pyim.io.fasta import read_fasta
from pyim.alignment.model import Alignment
from pyim.alignment.aligners.exact import ExactReadAligner
from pyim.alignment.aligners.inexact import WatermanAligner
from pyim.alignment.aligners.combined import CombinedReadAligner

logging.basicConfig(format='[%(asctime)s] %(levelname)-8s %(message)s', level=logging.INFO, datefmt='%H:%M:%S')


def main(fasta_file, vector_file, barcode_file, num_cores, read_limit):
    logging.info('Loading fasta sequences')
    fasta_seqs = _load_fasta(fasta_file, read_limit)

    # Align barcode sequences
    logging.info('Loading barcode definitions')
    barcode_seqs = list(read_fasta(barcode_file))

    aligner = CombinedReadAligner(exact_aligner=ExactReadAligner(), inexact_aligner=WatermanAligner())
    alignments = aligner.align_targets(fasta_seqs, barcode_seqs)

    # Load BLAT psl alignment output
    #pslReads = read_psl(pslInput)   # Just for reference
    #logging.info('Loading BLAT alignment')
    #bestReads = read_psl(args.read_alignment, onlyPutative=True)





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


def _parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--reads', required=True)
    parser.add_argument('--vectors', required=True)
    parser.add_argument('--barcodes', required=True)
    parser.add_argument('--read_alignment', required=True)
    parser.add_argument('--num_cores', type=int, default=1)
    parser.add_argument('--read_limit', type=int, default=None)
    return parser.parse_args()


if __name__ == '__main__':
    args = _parse_args()
    main(args.reads, args.vectors, args.barcodes, args.num_cores, args.read_limit)
