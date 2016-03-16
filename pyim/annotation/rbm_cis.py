from __future__ import absolute_import, division, print_function

#pylint: disable=wildcard-import,unused-wildcard-import,redefined-builtin
from builtins import *
#pylint: enable=wildcard-import,unused-wildcard-import,redefined-builtin

import logging
from os import path

import pandas as pd

#pylint: disable=import-error
from .rbm import rbm, WINDOW_PRESETS
#pylint: enable=import-error


def register(subparsers, name='rbm-cis'):
    parser = subparsers.add_parser(name, help=name + ' help')

    # Required arguments.
    parser.add_argument('input')
    parser.add_argument('output')
    parser.add_argument('--gtf', required=True)
    parser.add_argument('--cis_sites', required=True)

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--preset', choices=WINDOW_PRESETS.keys())
    group.add_argument('--window_sizes', nargs=4, type=int)

    # Optional arguments.
    parser.add_argument('--closest', default=False, action='store_true')

    # Set main for dispatch.
    parser.set_defaults(main=main)

    return parser


def main(args):
    logger = logging.getLogger()

    # Read insertions.
    insertions = pd.read_csv(args.input, sep='\t', dtype={'chrom': str})
    logger.info('Read %d insertions', len(insertions))

    # Read cis sites.
    cis_sites = pd.read_csv(args.cis_sites, sep='\t', dtype={'chrom': str})
    logger.info('Read %d cis sites', len(cis_sites))

    # Define windows.
    if args.window_sizes is not None:
        window_sizes = args.window_sizes
    else:
        window_sizes = WINDOW_PRESETS[args.preset]

    # Annotate cis sites.
    annotated_sites = rbm(cis_sites, args.gtf, window_sizes, logger,
                          closest=args.closest, verbose=True)

    # Extract and merge annotation with insertions.
    annotation = annotated_sites[['id', 'gene_id', 'gene_name',
                                  'gene_distance', 'window']]
    annotation = annotation.rename(columns={'id': 'cis_id'})

    annotated_ins = pd.merge(insertions, annotation, on='cis_id', how='left')

    # Write outputs.
    annotated_ins.to_csv(args.output, sep='\t', index=False)
    annotated_sites.to_csv(path.splitext(args.output)[0] + '.sites.txt',
                           sep='\t', index=False)
