from __future__ import absolute_import, division, print_function

#pylint: disable=wildcard-import,unused-wildcard-import,redefined-builtin
from builtins import *
#pylint: enable=wildcard-import,unused-wildcard-import,redefined-builtin

import itertools
import logging

import pandas as pd

from pyim.util.tabix import GtfFile

#pylint: disable=import-error
from ._model import Window
from ._util import select_closest
from .window import build_interval_trees, annotate_for_windows
#pylint: enable=import-error

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
    parser.add_argument('--closest', default=False, action='store_true')

    # Set main for dispatch.
    parser.set_defaults(main=main)

    return parser


def main(args):
    logger = logging.getLogger()

    # Read insertions.
    insertions = pd.read_csv(args.input, sep='\t', dtype={'chrom': str})
    logger.info('Read %d insertions', len(insertions))

    # Define windows.
    if args.window_sizes is not None:
        window_sizes = args.window_sizes
    else:
        window_sizes = WINDOW_PRESETS[args.preset]

    # Annotate insertions.
    annotated = rbm(insertions, args.gtf, window_sizes, logger,
                    closest=args.closest, verbose=True)
    annotated.to_csv(args.output, sep='\t', index=False)


def rbm(insertions, gtf_path, window_sizes, logger,
        closest=False, verbose=False):

    # Replace unstranded insertions with two stranded insertions.
    if (~insertions['strand'].isin({-1, 1})).any():
        logger.warning('Replacing unstranded insertions')
        converted = replace_unstranded(insertions)
    else:
        converted = insertions

    # Build annotation trees.
    logger.info('Building interval trees')
    gtf = GtfFile(gtf_path)
    trees = build_interval_trees(gtf)

    # Define windows.
    windows = build_windows(window_sizes)

    # Annotate insertions.
    logger.info('Annotating insertions')
    annotation = annotate_for_windows(
        converted, trees, windows, progress=verbose)

    if closest:
        logger.info('Reducing to closest features')
        annotation = select_closest(annotation, col='gene_distance')

    # Merge annotation with insertion frame.
    logger.info('Merging annotation')
    merged = pd.merge(insertions, annotation, on='id', how='left')

    return merged


def build_windows(window_sizes):
    us, ua, ds, da = window_sizes

    windows = [
        Window(0, 1, strand=1, incl_left=True, incl_right=True, name='is'),
        Window(0, 1, strand=-1, incl_left=True, incl_right=True, name='ia'),
        Window(-us, 0, strand=1, incl_left=True, incl_right=False, name='us'),
        Window(-ua, 0, strand=-1, incl_left=True, incl_right=False, name='ua'),
        Window(1, ds, strand=1, incl_left=False, incl_right=True, name='ds'),
        Window(1, da, strand=-1, incl_left=False, incl_right=True, name='da')]

    return windows


def replace_unstranded(insertions):
    """Replaces unstranded insertions with two stranded insertions."""

    # Split stranded and unstranded.
    mask = insertions['strand'].isin({-1, 1})
    stranded = insertions.ix[mask]
    unstranded = insertions.ix[~mask]

    # Convert unstranded into two stranded.
    converted = (_to_stranded(ins) for _, ins in unstranded.iterrows())
    converted = pd.DataFrame.from_records(
        itertools.chain.from_iterable(converted))

    return pd.concat((stranded, converted), ignore_index=True)


def _to_stranded(insertion):
    fwd = insertion.copy()
    fwd['strand'] = 1

    rev = insertion.copy()
    rev['strand'] = -1

    return [fwd, rev]
