import itertools

import numpy as np

from pyim.util.frozendict import frozendict


def select_closest(insertion, genes):
    """Filters potential hits for genes closest to the insertion."""

    if len(genes) == 0:
        return genes

    distances = np.array([_calc_distance(insertion, gene)
                          for gene in genes.itertuples()])  # yapf: disable
    abs_distances = np.abs(distances)

    return genes.loc[abs_distances == abs_distances.min()]


def filter_blacklist(genes, blacklist, field='gene_id'):
    """Filters potential hits against given blacklist."""
    return genes.loc[~genes[field].isin(blacklist)]


def annotate_insertion(insertion, hits):
    """Annotates insertion with given gene hits."""

    if len(hits) > 0:
        # Annotate insertion with overlapping genes.
        for row in hits.itertuples():
            gene_metadata = {
                'gene_id': row.gene_id,
                'gene_name': row.gene_name,
                'gene_distance': _calc_distance(insertion, row),
                'gene_orientation': 'sense' if row.strand \
                    == insertion.strand else 'antisense'
            }

            metadata = {**insertion.metadata,
                        **gene_metadata}  # yapf: disable
            yield insertion._replace(metadata=frozendict(metadata))
    else:
        # In case of no overlap, return original insertion.
        yield insertion


def _calc_distance(insertion, gene):
    assert not isinstance(gene.strand, str)

    if gene.start <= insertion.position < gene.end:
        return 0
    elif insertion.position > gene.end:
        dist = insertion.position - gene.end
    else:
        dist = insertion.position - gene.start

    dist *= gene.strand

    return dist
