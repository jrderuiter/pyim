from __future__ import absolute_import, division, print_function

#pylint: disable=wildcard-import,unused-wildcard-import,redefined-builtin
from builtins import *
#pylint: enable=wildcard-import,unused-wildcard-import,redefined-builtin

import itertools
import logging

import pandas as pd

#pylint: disable=import-error
from ..metadata import add_metadata
from ..filtering import filter_blacklist, select_closest
from .window import Window, annotate_windows
#pylint: enable=import-error

# Window format: (us, ua, ds, da)
WINDOW_PRESETS = {
    'SB': (20000, 10000, 25000, 5000),
    'MULV': (20000, 120000, 40000, 5000),
    'MMTV': (20000, 120000, 40000, 5000)
}


def annotate_rbm(insertions, gtf, window_preset=None, window_sizes=None):
    """Assigns insertions to genes using the rule-based-method (RBM) approach.

    Args:
        insertions (pandas.DataFrame): Insertions to annotate in DataFrame
            format. The frame is expected to contain at least the
            following columns: id, position, strand.
        gtf (str or GtfFile): Path to gtf file containing gene features.
            Alternatively, a GtfFile object may also be given instead of a path.
        window_preset (str): Preset to use for the RBM window sizes.
            Alternatively custom window sizes can be given using the
            *window_sizes* argument. Note that either *window_preset* or
            *window_sizes* must be provided.
        window_sizes (tuple[int]): Tuple of window sizes to use in the
            RBM mapping. Should specify four window sizes, for the following
            categories of insertions: upstream-sense, upstream-antisense,
            downstream-sense, downstream-antisense.

    Returns:
        pandas.DataFrame: Dataframe containing annotated insertions. Annotations
            are added as columns 'gene_id' and 'gene_name', which respectively contain the id and name of the annotated gene. An extra column
            'window' indicates which of the RBM windows was used for
            the annotation.

    """

    # Lookup windows.
    if window_preset is not None:
        window_sizes = WINDOW_PRESETS[window_preset]
    elif window_sizes is None:
        raise ValueError('Either window_sizes or window_preset must be given')

    # Replace unstranded insertions with two stranded insertions.
    if (~insertions['strand'].isin({-1, 1})).any():
        logging.warning('Replacing unstranded insertions')
        converted = _replace_unstranded(insertions)
    else:
        converted = insertions

    # Define windows.
    windows = _build_windows(window_sizes)

    # Annotate insertions.
    annotated = annotate_windows(converted, gtf, windows)

    return annotated


def _build_windows(window_sizes):
    us, ua, ds, da = window_sizes

    windows = [
        Window(0, 1, strand=1, incl_left=True, incl_right=True, name='is'),
        Window(0, 1, strand=-1, incl_left=True, incl_right=True, name='ia'),
        Window(-us, 0, strand=1, incl_left=True, incl_right=False, name='us'),
        Window(-ua, 0, strand=-1, incl_left=True, incl_right=False, name='ua'),
        Window(1, ds, strand=1, incl_left=False, incl_right=True, name='ds'),
        Window(1, da, strand=-1, incl_left=False, incl_right=True, name='da')]

    return windows


def _replace_unstranded(insertions):
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

    return (fwd, rev)


def register(subparsers, name='rbm'):
    """Registers the RBM annotator as a subparser."""

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
    parser.add_argument('--blacklist', default=None, nargs='+')

    # Set main for dispatch.
    parser.set_defaults(main=main)

    return parser


def main(args):
    """Main function for the RBM annotator command-line tool."""

    # Read insertions.
    insertions = pd.read_csv(args.input, sep='\t', dtype={'chrom': str})
    logging.info('Read %d insertions', len(insertions))

    # Annotate insertions.
    logging.info('Annotating insertions')
    annotated = annotate_rbm(insertions, args.gtf, window_preset=args.preset,
                             window_sizes=args.window_sizes)

    # Add metadata.
    logging.info('Adding annotation metadata')
    annotated = add_metadata(annotated, args.gtf)

    if args.blacklist is not None:
        logging.info('Filtering blacklisted genes')
        annotated = filter_blacklist(annotated, args.blacklist)

    if args.closest:
        logging.info('Selecting closest genes')
        annotated = select_closest(annotated)


    annotated.to_csv(args.output, sep='\t', index=False)
