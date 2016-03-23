from __future__ import absolute_import, division, print_function

#pylint: disable=wildcard-import,unused-wildcard-import,redefined-builtin
from builtins import *
#pylint: enable=wildcard-import,unused-wildcard-import,redefined-builtin

import itertools
import logging

import pandas as pd
from intervaltree import IntervalTree
from tqdm import tqdm

from pyim.util.tabix import GtfFile

# pylint: disable=import-error
from ..filtering import filter_blacklist, select_closest
from ..util import build_interval_trees, numeric_strand
# pylint: enable=import-error

def annotate_windows(insertions, gtf, windows):
    """Assigns insertions to genes that fall within the given windows.

    Args:
        insertions (pandas.DataFrame): Insertions to annotate in DataFrame
            format. The frame is expected to contain at least the
            following columns: id, position, strand.
        gtf (str or GtfFile): Path to gtf file containing gene features.
            Alternatively, a GtfFile object may also be given instead of a path.
        windows (list[Window]): List of windows to inspect for genes.

    Returns:
        pandas.DataFrame: Dataframe containing annotated insertions. Annotations
            are added as columns 'gene_id' and 'gene_name', which respectively contain the id and name of the annotated gene. An extra column
            'window' indicates which of the RBM windows was used for
            the annotation.

    """

    if isinstance(gtf, str):
        gtf = GtfFile(gtf)

    # Build lookup trees.
    trees = build_interval_trees(gtf)

    # Generate queries (insertion/window combinations).
    ins_gen = (row for _, row in insertions.iterrows())
    queries = itertools.product(ins_gen, windows)

    queries = tqdm(queries, unit='query',
                   total=len(insertions) * len(windows))

    # Generate annotation for queries and merge into frame.
    annotations = (_annotate_window(ins, window, trees)
                   for ins, window in queries)
    annotation = pd.concat(annotations, ignore_index=True)

    # Merge annotation with insertions.
    annotated = pd.merge(insertions, annotation, on='id', how='left')

    return annotated


class Window(object):
    """Class representing a (relative) window to inspect for genes.

    The window may be an actual window corresponding to a real chromosome
    location, in which case start and end represent the actual window
    boundaries, and reference and strand represent the actual chromosome
    and strand of the window.

    Alternatively, the window may also represent a relative window. In this
    case start is typically negative and end is typically positive, whilst
    reference is typically omitted and strand is optional. This relative window
    can be applied to an actual position using the apply method, which
    effectively calculates the given window around that position.

    Args:
        start (int): Start of the window.
        end (int): End of the window.
        reference (str): Chromosome of the window (optional).
        strand (int): Relative strand of window (optional).
        incl_left (bool): Whether to include partially (left)
            overlapping features.
        incl_right (bool): Whether to include partially (right)
            overlapping features.

    """

    def __init__(self, start, end, reference=None, strand=None,
                 incl_left=True, incl_right=True, name=None):
        self.reference = reference
        self.start = start
        self.end = end
        self.strand = strand

        self.incl_left = incl_left
        self.incl_right = incl_right

        self.name = name

    def apply(self, reference, location, strand):
        """Applies a relative window to specific location and strand.

        For example, a relative window of Window(start=-1000, end=1000,
        strand=-1) applied to position (2, 3000, -1) will become
        Window(ref=2, start=2000, end=4000, strand=1).

        Args:
            reference (str): Chromosome name of the reference position.
            location (int): Reference genomic position.
            strand (int): Reference genomic strand.

        """

        # Determine start/end position.
        if strand == 1:
            start = location + self.start
            end = location + self.end

            incl_left = self.incl_left
            incl_right = self.incl_right
        elif strand == -1:
            start = location - self.end
            end = location - self.start

            incl_right = self.incl_left
            incl_left = self.incl_right
        else:
            raise ValueError('Unknown value for strand ({})'
                             .format(strand))

        # Determine new strand.
        if self.strand is not None:
            new_strand = self.strand * strand
        else:
            new_strand = None

        return Window(start, end, reference, new_strand,
                      incl_left, incl_right, name=self.name)


def _annotate_window(insertion, window, feature_trees):
    """Annotates insertion for features in trees using given window."""

    # Apply window for insertion.
    applied_window = window.apply(
        insertion['chrom'], insertion['position'], insertion['strand'])

    # Fetch features within window.
    features = _fetch_in_window(feature_trees, applied_window)

    # Extract feature values.
    frame = pd.DataFrame.from_records(
        ({'id': insertion['id'],
          'gene_id': feature['gene_id'],
          'gene_name': feature['gene_name']}
         for feature in features))

    # Include window name if known.
    if window.name is not None:
        frame['window'] = window.name

    return frame


def _fetch_in_window(trees, window):
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
    # Read annotation.
    insertions = pd.read_csv(args.input, sep='\t', dtype={'chrom': str})
    logging.info('Read %d insertions', len(insertions))

    # Define windows.
    logging.info('Annotating insertions')
    half_size = args.window_size // 2
    window = Window(start=-half_size, end=half_size)

    # Annotate insertions.
    annotated = annotate_windows(insertions, args.gtf, [window])

    if args.blacklist is not None:
        logging.info('Filtering blacklisted genes')
        annotated = filter_blacklist(annotated, args.blacklist)

    if args.closest:
        logging.info('Selecting closest genes')
        annotated = select_closest(annotated)

    # Merge annotation.
    annotated.to_csv(args.output, sep='\t', index=False)
