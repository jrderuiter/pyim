from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
from builtins import (ascii, bytes, chr, dict, filter, hex, input,
                      int, map, next, oct, open, pow, range, round,
                      str, super, zip)
from future.utils import native_str

from argparse import ArgumentParser
from os import path

import logging
import pandas as pd

from pyim.cis.cimpl import map_insertions
from pyim.util.insertions import subset_samples

from ._logging import print_header, print_footer


def setup_parser():
    parser = ArgumentParser(prog='pyim-cis')

    parser.add_argument('input')
    parser.add_argument('output')

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--pattern', default=None)
    group.add_argument('--system', choices={'SB'}, default=None)

    parser.add_argument('--genome', choices={'mm10'}, default='mm10')
    parser.add_argument('--chromosomes', nargs='+', default=None)
    parser.add_argument('--scales', nargs='+', type=int, default=30000)
    parser.add_argument('--samples', nargs='+', default=None)

    parser.add_argument('--iterations', type=int, default=1000)
    parser.add_argument('--lhc_method', choices={'none', 'exclude'},
                        default='exclude')

    # parser.add_argument('--strand_homogeneity', type=float, default=0.75)

    parser.add_argument('--alpha', type=float, default=0.05)

    parser.add_argument('--threads', type=int, default=1)
    parser.add_argument('--verbose', default=False, action='store_true')

    return parser


def main():
    logger = logging.getLogger()

    # Parse arguments.
    parser = setup_parser()
    args = parser.parse_args()

    # Print header.
    print_header(logger, command='cis')

    # Read insertions.
    insertions = pd.read_csv(args.input, sep=native_str('\t'),
                             dtype={'chrom': str})
    logger.info('Read {} insertions'.format(len(insertions)))

    # Subset to samples if needed.
    if args.samples is not None:
        logger.info('Subsetting to {} samples'.format(len(args.samples)))
        insertions = subset_samples(insertions, args.samples, logger=logger)

    # Run cimpl on insertions.
    logger.info('Running CIMPL in R')

    cis, mapping = map_insertions(
        insertions, scales=args.scales, genome=args.genome, alpha=args.alpha,
        system=args.system, pattern=args.pattern, lhc_method=args.lhc_method,
        chromosomes=args.chromosomes, iterations=args.iterations,
        threads=args.threads, verbose=args.verbose)

    # Annotate insertions with cis mapping.
    logger.info('Merging CIMPL annotation')

    mapping_tmp = mapping.rename(columns={'insertion_id': 'id'})
    insertions = pd.merge(insertions, mapping_tmp, on='id')

    # Write out outputs.
    logger.info('Writing outputs')

    cis_path = path.splitext(args.output)[0] + '.sites.txt'
    cis.to_csv(cis_path, sep=native_str('\t'), index=False)

    insertions.to_csv(args.output, sep=native_str('\t'), index=False)

    print_footer(logger)


if __name__ == '__main__':
    main()
