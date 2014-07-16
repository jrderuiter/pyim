#! /usr/bin/env python

import argparse
from os import path
import logging

from pyim.io import makedirs_safe
from pyim.pipeline.sb.pipeline import sb_pipeline
from pyim.alignment.vector.aligners import ExactReadAligner, TruncatedTargetAligner, ExonerateReadAligner, \
                                           ChainedReadAligner, ParallelTargetAligner, ParallelReadAligner
from pyim.alignment.vector.filters import QueryEndFilter, MismatchFilter, BestScoreFilter


DATA_DIR = path.normpath(path.join(path.dirname(__file__), '../../data'))
DEFAULT_CONTAMINANTS = path.join(DATA_DIR, 'SB_contaminants.fa')
DEFAULT_BARCODES     = path.join(DATA_DIR, 'SB_barcodes.fa')
DEFAULT_VECTORS      = path.join(DATA_DIR, 'vec.fa')


def main(options):

    # Create logger
    logger = logging.getLogger('PyIM')
    logger.setLevel(logging.DEBUG)

    ch = logging.StreamHandler()
    ch.setLevel(logging.DEBUG)

    formatter = logging.Formatter('%(asctime)s - %(name)s - %(message)s')
    ch.setFormatter(formatter)
    logger.addHandler(ch)

    # Create output dir
    makedirs_safe(options.output_dir)

    # Setup aligners
    sb_aligner, t7_aligner, bc_aligner = _setup_aligners(options.threads)
    vector_aligners = {'sb': sb_aligner, 't7': t7_aligner, 'bc': bc_aligner}

    # Run SB pipeline
    insertions = sb_pipeline(reads_file=options.reads_file, run_name=options.run_name,
                             reference=options.reference, vector_aligners=vector_aligners,
                             vector_file=options.vector_file, barcode_file=options.barcode_file,
                             barcode_sample_file=options.barcode_sample_file, work_dir=options.output_dir,
                             contaminant_file=options.contaminant_file, threads=options.threads,
                             min_sequence_len=options.min_seq_len, min_ulp=options.min_ulp)

    # Write output
    output_file = path.join(options.output_dir, 'insertions.txt')
    insertions.to_csv(output_file, sep='\t', index=False)


def _parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument('-i', '--input_reads',      dest='reads_file', required=True)
    parser.add_argument('-o', '--output_dir',       required=True)

    parser.add_argument('-r', '--reference',        required=True)
    parser.add_argument('-m', '--barcode-mapping',  dest='barcode_sample_file', required=True)
    parser.add_argument('-n', '--run-name',         dest='run_name', required=True)

    parser.add_argument('-v', '--vector-file',      dest='vector_file',      default=DEFAULT_VECTORS)
    parser.add_argument('-b', '--barcode-file',     dest='barcode_file',     default=DEFAULT_BARCODES)
    parser.add_argument('-c', '--contaminant-file', dest='contaminant_file', default=DEFAULT_CONTAMINANTS)

    parser.add_argument('-l', '--min-seq-len',      dest='min_seq_len', default=15, type=int)
    parser.add_argument('-u', '--min-ulp',          dest='min_ulp',     default=2, type=int)

    parser.add_argument('-t', '--threads',          default=1, type=int)

    return parser.parse_args()


def _setup_aligners(threads):
    trunc1_aligner = TruncatedTargetAligner(ExactReadAligner(), 1)
    trunc3_aligner = TruncatedTargetAligner(ExactReadAligner(), 3)
    ex_aligner = ExonerateReadAligner(min_score=60, filters=[QueryEndFilter(verbose=True)])

    exhaust_t7_filters = [QueryEndFilter(verbose=True),
                          MismatchFilter(2, include_start=True, verbose=True),
                          BestScoreFilter(verbose=True)]
    exhaust_aligner = ExonerateReadAligner(min_score=40, exhaustive=True, bestn=10)
    exhaust_aligner = ParallelReadAligner(exhaust_aligner, threads, filters=exhaust_t7_filters)

    t7_aligner = ChainedReadAligner([ExactReadAligner(), trunc1_aligner, trunc3_aligner, ex_aligner, exhaust_aligner])
    bc_aligner = ParallelTargetAligner(ExactReadAligner(), threads)
    sb_aligner = ChainedReadAligner([ExactReadAligner(), ExonerateReadAligner(min_score=60)])

    return sb_aligner, t7_aligner, bc_aligner


if __name__ == '__main__':
    main(_parse_args())
