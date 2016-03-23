import itertools
from intervaltree import IntervalTree


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


def numeric_strand(strand):
    """Converts strand to its numeric (integer) representation."""

    if isinstance(strand, int):
        return strand
    elif isinstance(strand, float):
        return int(strand)
    else:
        if strand == '+':
            return 1
        elif strand == '-':
            return -1
        else:
            raise ValueError('Unknown value {} for strand'.format(strand))
