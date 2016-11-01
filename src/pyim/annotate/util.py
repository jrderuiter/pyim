import itertools
from pathlib import Path

from intervaltree import IntervalTree
import numpy as np

from pyim.util.tabix import GtfFile


def build_interval_trees(reference_gtf):
    """Builds an interval tree of genes for each chromosome in gtf."""

    if isinstance(reference_gtf, (str, Path)):
        reference_gtf = GtfFile(reference_gtf)

    # Only select gene features for now.
    genes = reference_gtf.fetch(filters={'feature': 'gene'})

    # Note, below code assumes that genes are ordered by contig.

    trees = {}
    for contig, grp in itertools.groupby(genes, lambda r: r['contig']):
        # Build a tree for each individual chromosome.
        intervals = ((g['start'], g['end'], dict(g)) for g in grp
                     if g['end'] > g['start'])  # Avoid null intervals.
        trees[contig] = IntervalTree.from_tuples(intervals)

    return trees


def numeric_strand(strand):
    """Converts strand to its numeric (integer) representation."""

    if isinstance(strand, int):
        return strand
    elif isinstance(strand, (float, np.generic)):
        return int(strand)
    elif isinstance(strand, str):
        if strand == '+':
            return 1
        elif strand == '-':
            return -1

    raise ValueError('Unknown value {} for strand (type: {})'
                     .format(strand, type(strand)))
