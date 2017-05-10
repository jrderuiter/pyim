import itertools
import operator
import sys

import numpy as np


def select_closest(insertions, field='gene_distance'):
    """Selects genes that are closest to the annotated insertions.

    Parameters:
        insertions (iterable[Insertion]): Annotated insertions that are to
            be filtered. The frame is expected to contain at least the
            following columns: id, position, strand, *dist_col*.
        field (str): Name of the column containing the distance to
            the gene or feature. Can be added using the add_metadata function.

    Returns:
        iterable[Insertion]: Filtered annotated insertions, which have been
            reduced to only include the genes closest to the insertions.

    """

    # Group insertions by id.
    id_getter = operator.attrgetter('id')
    insertions = sorted(insertions, key=id_getter)
    grouped = itertools.groupby(insertions, key=id_getter)

    for _, group in grouped:
        group = list(group)

        # Yield closest insertions (with minimum distance).
        dists = np.abs([ins.metadata.get(field, sys.maxsize) for ins in group])
        yield from itertools.compress(group, dists == dists.min())


def filter_blacklist(insertions, blacklist, field='gene_name'):
    """Filters annotations that assign insertions to blacklisted genes.

    Args:
        insertions (iterable[Insertion]):
        blacklist (iterable[str]): List of blacklisted gene ids to filter.
        field (str): Name of the column containing the id of the genes.

    Returns:
        iterable[Insertion]: Filtered annotated insertions, which have been
            reduced to remove blacklisted genes.

    """

    # Ensure blacklist is a set.
    blacklist = set(blacklist)

    # Drop any genes with a gene id in the blacklist.
    for insertion in insertions:
        if not insertion.metadata.get(field, None) in blacklist:
            yield insertion
