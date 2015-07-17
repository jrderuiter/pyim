from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
from builtins import (ascii, bytes, chr, dict, filter, hex, input,
                      int, map, next, oct, open, pow, range, round,
                      str, super, zip)
from future.utils import native_str

from argparse import ArgumentParser
from pathlib import Path

import numpy as np
import pandas as pd
from toolz import curry

from pyim.cis.cimpl import cimpl, get_cis, get_cis_mapping


def setup_parser():
    parser = ArgumentParser(prog='pyim-cis')

    parser.add_argument('input', type=Path)
    parser.add_argument('output', type=Path)

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--pattern', default=None)
    group.add_argument('--system', choices={'SB'}, default=None)

    parser.add_argument('--genome', choices={'mm10'}, default='mm10')
    parser.add_argument('--chromosomes', nargs='+', default=None)
    parser.add_argument('--scales', nargs='+', type=int, default=30000)

    parser.add_argument('--strand_homogeneity', type=float, default=0.75)

    parser.add_argument('--alpha', type=float, default=0.05)
    parser.add_argument('--iterations', type=int, default=1000)
    parser.add_argument('--lhc_method', choices={'none', 'exclude'},
                        default='exclude')

    parser.add_argument('--threads', type=int, default=1)
    parser.add_argument('--verbose', default=False, action='store_true')

    return parser


def main():
    parser = setup_parser()
    args = parser.parse_args()

    # Read frame.
    ins_frame = pd.read_csv(str(args.input), sep=native_str('\t'))

    # Run cimpl.
    cimpl_obj = cimpl(ins_frame, scales=args.scales, genome=args.genome,
                      system=args.system, pattern=args.pattern,
                      lhc_method=args.lhc_method, chromosomes=args.chromosomes,
                      iterations=args.iterations, threads=args.threads,
                      verbose=args.verbose)

    # Extract cis and cis mapping from object.
    cis = get_cis(cimpl_obj, alpha=args.alpha, mul_test=True)
    cis_mapping = get_cis_mapping(cimpl_obj, cis_frame=cis)

    # Annotate insertions with cis mapping.
    ins_annotated = pd.merge(ins_frame, cis_mapping, on='insertion_id')

    # Determine strand of cis sites.
    strand_func = curry(_strandedness, min_homogeneity=args.strand_homogeneity)
    cis_strand = ins_annotated.groupby('cis_id').apply(strand_func)

    # Merge strand information with cis sites.
    cis = pd.merge(cis, cis_strand.reset_index(), on='cis_id')

    # Rename and reshuffle cis columns.
    cis = cis.rename(columns={'peak_location': 'location',
                              'peak_height': 'height'})
    cis = cis[['cis_id', 'seqname', 'location', 'strand', 'scale',
               'n_insertions', 'p_value', 'start', 'end', 'height', 'width',
               'strand_mean', 'strand_homogeneity']]

    # Write out outputs.
    cis.to_csv(str(args.output.with_suffix('.sites.txt')),
               sep=native_str('\t'), index=False)
    ins_annotated.to_csv(str(args.output), sep=native_str('\t'), index=False)


def _strandedness(insertions, min_homogeneity):
    strand_mean = insertions.strand.mean()
    strand = int(np.sign(strand_mean))

    if strand != 0:
        homogeneity = (insertions.strand == strand).sum() / len(insertions)
    else:
        homogeneity = 0.5

    if homogeneity < min_homogeneity:
        strand = 0

    return pd.Series(dict(strand=strand,
                          strand_mean=strand_mean,
                          strand_homogeneity=homogeneity))



if __name__ == '__main__':
    main()
