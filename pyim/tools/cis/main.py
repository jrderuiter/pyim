import argparse

import numpy as np
import pandas as pd

from pyim.tools.cis import cimpl


def cis_main(args):
    insertions = pd.read_csv(args.input, sep='\t')

    # Run cimpl on our insertions.
    cimpl_obj = cimpl.cimpl(insertions, scales=args.scales, n_iterations=args.n_iter,
                            chromosomes=args.chromosomes, system=args.system,
                            specificity_pattern=args.specificity_pattern,
                            genome=args.genome, threads=args.threads, verbose=True)

    # Extract required information from the cimpl object.
    cis = cimpl.cis(cimpl_obj, alpha=args.alpha, mul_test=True)
    cis_mapping = cimpl.cis_mapping(cimpl_obj, cis)

    # Add strand + homogeneity to cis frame.
    cis_strand = cis_strandedness(insertions, cis_mapping)
    cis = pd.merge(cis, cis_strand, on='id')

    # Rename some of the columns for conformity.
    cis.rename(columns={'chromosome': 'seqname',
                        'peak_location': 'location',
                        'peak_height': 'height'}, inplace=True)

    # Write output files!
    cis.to_csv(args.output_base + '.cis.txt', sep='\t', index=False)
    cis_mapping.to_csv(args.output_base + '.cis.mapping.txt', index=False, sep='\t')


def cis_strandedness(insertions, cis_mapping):
    cis_merged = pd.merge(insertions, cis_mapping, left_on='name', right_on='insertion_id')
    strand_mean = cis_merged.groupby('cis_id')['strand'].mean()

    # Decide on strand, replace zero's with 1 (fwd strand)..
    cis_strand = np.sign(strand_mean).astype(int)
    cis_strand.ix[cis_strand == 0] = 1

    # Calculate the homogeneity of the cis.
    cis_homogeneity = (strand_mean + 1) / 2
    cis_homogeneity = np.maximum(cis_homogeneity, 1 - cis_homogeneity)
    cis_homogeneity.name = 'strand_homogeneity'

    # Merge results!
    frame = pd.concat([cis_strand, cis_homogeneity], axis=1).reset_index()
    frame = frame.rename(columns={'cis_id': 'id'})

    return frame


def _parse_args():
    # TODO: allow self-supplied reference?

    parser = argparse.ArgumentParser()

    parser.add_argument('input')
    parser.add_argument('output_base')

    parser.add_argument('--genome', default='mm10', choices=['mm10'])
    parser.add_argument('--scales', default=30000, type=int, nargs='+')
    parser.add_argument('--n-iter', default=1000, type=int)
    parser.add_argument('--alpha', default=0.05, type=float)
    parser.add_argument('--chromosomes', default=None, nargs='+')
    parser.add_argument('--threads', default=1, type=int)

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--system', choices=['SB', 'PB', 'MMTV', 'MuLV'])
    group.add_argument('--specificity-pattern')

    return parser.parse_args()


def main():
    cis_main(_parse_args())


if __name__ == '__main__':
    main()
