from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
from builtins import (ascii, bytes, chr, dict, hex, input,
                      int, map, next, oct, open, pow, range, round,
                      str, super, zip)  # filter

import logging

import pandas as pd

from pyim.util.tabix import GtfFile

from ._model import Window
from ._util import select_closest
from .window import build_interval_trees, annotate_for_windows


# Window format: (us, ua, ds, da)
WINDOW_PRESETS = {
    'SB': (20000, 10000, 25000, 5000),
    'MULV': (20000, 120000, 40000, 5000),
    'MMTV': (20000, 120000, 40000, 5000)
}


def register(subparsers, name='rbm'):
    parser = subparsers.add_parser(name, help=name + ' help')

    # Required arguments.
    parser.add_argument('input')
    parser.add_argument('output')
    parser.add_argument('--gtf', required=True)

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--preset', choices=WINDOW_PRESETS.keys())
    group.add_argument('--window_sizes', nargs=4, type=int)

    # Optional arguments.
    # parser.add_argument('--feature_type', default='gene',
    #                     choices={'gene', 'transcript'})
    # parser.add_argument('--id_column', default='insertion_id')
    parser.add_argument('--closest', default=False, action='store_true')

    # Set main for dispatch.
    parser.set_defaults(main=main)

    return parser


def main(args):
    logger = logging.getLogger()

    # Read insertions.
    insertions = pd.read_csv(args.input, sep='\t', dtype={'chrom': str})
    logger.info('Read {} insertions'.format(len(insertions)))

    # Build annotation trees.
    logger.info('Building interval trees')
    gtf = GtfFile(args.gtf)
    trees = build_interval_trees(gtf)

    # Define windows.
    if args.preset is not None:
        window_sizes = WINDOW_PRESETS[args.preset]
    else:
        window_sizes = args.window_sizes

    windows = build_windows(window_sizes)

    # Annotate insertions.
    logger.info('Annotating insertions')
    annotation = annotate_for_windows(
        insertions, trees, windows, progress=True)

    if args.closest:
        logger.info('Reducing to closest features')
        annotation = select_closest(annotation, col='gene_distance')

    # Merge annotation with insertion frame.
    logger.info('Merging annotation')
    merged = pd.merge(insertions, annotation, on='id', how='left')
    merged.to_csv(args.output, sep='\t', index=False)


def build_windows(ranges):
    us, ua, ds, da = ranges

    windows = [
        Window(0, 1, strand=1, incl_left=True, incl_right=True, name='is'),
        Window(0, 1, strand=-1, incl_left=True, incl_right=True, name='ia'),
        Window(-us, 0, strand=1, incl_left=True, incl_right=False, name='us'),
        Window(-ua, 0, strand=-1, incl_left=True, incl_right=False, name='ua'),
        Window(1, ds, strand=1, incl_left=False, incl_right=True, name='ds'),
        Window(1, da, strand=-1, incl_left=False, incl_right=True, name='da')
    ]

    return windows
