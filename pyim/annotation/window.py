from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
from builtins import (ascii, bytes, chr, dict, filter, hex, input,
                      int, map, next, oct, open, pow, range, round,
                      str, super, zip)

import itertools
import logging

import pandas as pd
from intervaltree import IntervalTree
from tqdm import tqdm

from pyim.util.tabix import GtfFile

from ._model import Window
from ._util import feature_distance, numeric_strand, select_closest


def register(subparsers, name='window'):
    parser = subparsers.add_parser(name, help=name + ' help')

    # Required arguments.
    parser.add_argument('input')
    parser.add_argument('output')
    parser.add_argument('--gtf', required=True)

    # Optional arguments.
    parser.add_argument('--closest', default=False, action='store_true')
    parser.add_argument('--window_size', default=20000, type=int)

    # Set main for dispatch.
    parser.set_defaults(main=main)

    return parser


def main(args):
    logger = logging.getLogger()

    # Read annotation.
    insertions = pd.read_csv(args.input, sep='\t', dtype={'chrom': str})
    logger.info('Read {} insertions'.format(len(insertions)))

    # Build lookup trees.
    logger.info('Building interval trees')
    gtf = GtfFile(args.gtf)
    trees = build_interval_trees(gtf)

    # Define windows.
    logger.info('Annotating insertions')
    half_size = args.window_size // 2
    window = Window(start=-half_size, end=half_size)

    # Annotate insertions.
    annotation = annotate_for_windows(
        insertions, trees, [window], progress=True)

    if args.closest:
        # Sub-select for closest features.
        logger.info('Reducing to closest features')
        annotation = select_closest(annotation, col='gene_distance')

    # Merge annotation.
    logger.info('Merging annotation')
    merged = pd.merge(insertions, annotation, on='id', how='left')
    merged.to_csv(args.output, sep='\t', index=False)


def annotate_for_windows(insertions, trees, windows, progress=False):
    """Annotates insertions for features in trees using given windows."""

    # Generate queries (insertion/window combinations).
    ins_gen = (row for _, row in insertions.iterrows())
    queries = itertools.product(ins_gen, windows)

    if progress:
        queries = tqdm(queries, unit='query',
                       total=len(insertions) * len(windows))

    # Generate annotation for queries.
    annotations = (_annotate_for_window(ins, trees, window)
                   for ins, window in queries)

    # Merge annotations into single frame.
    annotation = pd.concat(annotations, ignore_index=True)

    return annotation


def _annotate_for_window(insertion, trees, window):
    """Annotates insertion for features in trees using given window."""

    # Apply window for insertion.
    applied_window = window.apply(
        insertion['chrom'], insertion['position'], insertion['strand'])

    # Fetch features within window.
    features = fetch_in_window(trees, applied_window)

    # Extract feature values.
    values = ((f['gene_id'],
               f['gene_name'],
               feature_distance(f, insertion['position']))
              for f in features)

    try:
        id_, name, distance = zip(*values)
    except ValueError:
        id_, name, distance = [], [], []

    # Convert to frame.
    frame = pd.DataFrame({
        'id': insertion['id'],
        'gene_id': id_,
        'gene_name': name,
        'gene_distance': distance})

    # Include window name if known.
    if window.name is not None:
        frame['window'] = window.name

    return frame


def fetch_in_window(trees, window):
    """Fetches features within given window in the interval trees."""

    # Find overlapping features.
    try:
        tree = trees[window.reference]
        overlap = tree[window.start:window.end]
    except KeyError:
        overlap = []

    # Extract features.
    features = (interval[2] for interval in overlap)

    # Filter inclusive/exclusive if needed.
    if not window.incl_left:
        features = (f for f in features if f['start'] > window.start)

    if not window.incl_right:
        features = (f for f in features if f['end'] < window.end)

    # Filter for strand if needed.
    if window.strand is not None:
        features = (f for f in features
                    if numeric_strand(f['strand']) == window.strand)

    return list(features)


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
