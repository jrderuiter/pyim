#! /usr/bin/env python

import argparse
from os import path
import logging


from pyim.io import makedirs_safe
from pyim.pipeline.sb.pipeline import sb_pipeline
from pyim.alignment.vector.aligners import ExactReadAligner, TruncatedTargetAligner, ExonerateReadAligner, \
                                           ChainedReadAligner, ParallelTargetAligner, ParallelReadAligner
from pyim.alignment.vector.filters import QueryEndFilter, MismatchFilter, BestScoreFilter


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


def _parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument('-i', '--input_reads', dest='reads_file', required=True)
    parser.add_argument('-o', '--output_file', dest='output_file', required=True)

    parser.add_argument('-r', '--reference', dest='reference', required=True)
    parser.add_argument('-v', '--vector_file', dest='vector_file', required=True)
    parser.add_argument('-b', '--barcode_file', dest='barcode_file', required=True)
    parser.add_argument('-m', '--barcode_mapping', dest='barcode_sample_file', required=True)
    parser.add_argument('-c', '--contaminant_file', dest='contaminant_file', default=None)

    parser.add_argument('-l', '--min-seq-len', dest='min_seq_len', default=15, type=int)
    parser.add_argument('-u', '--min-ulp', dest='min_ulp', default=2, type=int)

    parser.add_argument('-t', '--threads', dest='threads', default=1, type=int)

    return parser.parse_args()


def main(options):

    # create logger
    logger = logging.getLogger('PyIM')
    logger.setLevel(logging.DEBUG)

    ch = logging.StreamHandler()
    ch.setLevel(logging.DEBUG)

    formatter = logging.Formatter('%(asctime)s - %(name)s - %(message)s')
    ch.setFormatter(formatter)
    logger.addHandler(ch)

    # Create output dir
    output_dir = path.dirname(options.output_file)
    makedirs_safe(output_dir)

    # Setup aligners
    sb_aligner, t7_aligner, bc_aligner = _setup_aligners(options.threads)
    vector_aligners = {'sb': sb_aligner, 't7': t7_aligner, 'bc': bc_aligner}

    # Run SB pipeline
    insertions = sb_pipeline(options.reads_file, options.reference, output_dir, vector_aligners,
                             options.vector_file, options.barcode_file,
                             options.barcode_sample_file, options.contaminant_file,
                             options.threads, options.min_seq_len, options.min_ulp)

    # Write output
    insertions.to_csv(options.output_file, sep='\t', index=False)


if __name__ == '__main__':
    main(_parse_args())

