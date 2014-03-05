#! /usr/bin/env python

import os
import argparse

from pyim.sb.pipeline import sb_pipeline
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


if __name__ == '__main__':
    args = _parse_args()

    sb_aligner, t7_aligner, bc_aligner = _setup_aligners(args.threads)

    insertions = sb_pipeline(args.reads_file, args.reference, args.vector_file, args.barcode_file,
                             args.barcode_sample_file, args.contaminant_file, args.threads,
                             args.min_seq_len, args.min_ulp, sb_aligner=sb_aligner,
                             t7_aligner=t7_aligner, bc_aligner=bc_aligner)

    output_dir = os.path.dirname(args.output_file)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)

    insertions.to_csv(args.output_file, sep='\t', index=False)
