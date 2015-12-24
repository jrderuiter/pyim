from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
from builtins import (ascii, bytes, chr, dict, filter, hex, input,
                      int, map, next, oct, open, pow, range, round,
                      str, super, zip)

import itertools
import logging

import pandas as pd
from intervaltree import IntervalTree

from pyim.util.tabix import GtfFile

from ._model import Window


def register(subparsers, name='window'):
    parser = subparsers.add_parser(name, help=name + ' help')

    # Required arguments.
    parser.add_argument('input')
    parser.add_argument('output')
    parser.add_argument('--gtf', required=True)

    # Optional arguments.
    # parser.add_argument('--feature_type', default='gene')
    parser.add_argument('--window_size', default=20000, type=int)

    # Set main for dispatch.
    parser.set_defaults(main=main)

    return parser


def main(args):
    logger = logging.getLogger()

    insertions = pd.read_csv(args.input, sep='\t', dtype={'chrom': str})
    logger.info('Read {} insertions'.format(len(insertions)))

    logger.info('Building interval trees')
    gtf = GtfFile(args.gtf)
    trees = build_interval_trees(gtf)

    logger.info('Annotating insertions')
    half_size = args.window_size // 2
    window = Window(start=-half_size, end=half_size)

    annotation = annotate_for_windows(insertions, trees, [window])

    logger.info('Merging annotation')
    merged = pd.merge(insertions, annotation, on='id', how='left')
    merged.to_csv(args.output, sep='\t', index=False)


def annotate_for_windows(insertions, trees, windows):
    """Annotates insertions for features in trees using given windows."""

    if isinstance(insertions, pd.DataFrame):
        insertions = (row for _, row in insertions.iterrows())

    queries = itertools.product(insertions, windows)

    annotation = pd.concat((_annotate_for_window(ins, trees, window)
                           for ins, window in queries), ignore_index=True)

    return annotation


def _annotate_for_window(insertion, trees, window):
    """Annotates insertion for features in trees using given window."""

    # Apply window for insertion.
    applied_window = window.apply(
        insertion['chrom'], insertion['position'], insertion['strand'])

    # Fetch features within window.
    features = fetch_in_window(trees, applied_window)

    # Convert to frame.
    frame = pd.DataFrame({
        'id': insertion['id'],
        'gene_name': [f['gene_name'] for f in features]})

    # Include window name if known.
    if window.name is not None:
        frame['window'] = window.name

    return frame


def fetch_in_window(trees, window):
    """Fetches features within given window in the interval trees."""

    if window.strand is not None:
        raise NotImplementedError()

    try:
        tree = trees[window.reference]
        overlap = tree[window.start:window.end]
    except KeyError:
        overlap = []

    features = [interval[2] for interval in overlap]

    if window.strand is not None:
        features = [f for f in features
                    if _strand_numeric(f['strand']) == window.strand]

    return features


def _strand_numeric(strand):
    return 1 if strand == '+' else -1


def build_interval_trees(gtf):
    """Builds an interval tree of genes for each chromosome in gtf."""

    # Only select gene features for now.
    genes = gtf.fetch(filters={'feature': 'gene'})

    trees = {}
    for contig, grp in itertools.groupby(genes, lambda r: r.contig):
        # Build a tree for each individual chromosome.
        intervals = ((g.start, g.end, dict(g)) for g in grp
                      if g.end > g.start)  # Avoid null intervals.
        trees[contig] = IntervalTree.from_tuples(intervals)

    return trees
