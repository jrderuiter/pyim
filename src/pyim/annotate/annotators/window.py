from collections import namedtuple

import pandas as pd

from .base import Annotator


class WindowAnnotator(Annotator):
    """Window annotator class."""

    def __init__(self, windows):
        super().__init__()
        self._windows = windows

    @classmethod
    def from_window_size(cls, window_size, **kwargs):
        """Creates instance using given window size."""

        window = Window(
            upstream=window_size,
            downstream=window_size,
            strand=None,
            name='default',
            strict_left=False,
            strict_right=False)

        return cls(windows=[window], **kwargs)

    def match(self, insertions, genes):
        def _match_gen(insertions):
            for row in insertions.itertuples():
                for window in self._windows:
                    region = window.apply(row.chromosome, row.position,
                                          row.strand)
                    gene_ids = self._get_overlap(genes, region)

                    for gene_id in gene_ids:
                        yield (row.id, gene_id, window.name)

        return pd.DataFrame.from_records(
            _match_gen(insertions), columns=['id', 'gene_id', 'window'])

    def _get_overlap(self, genes, region):
        try:
            overlap = genes.gloc.search(
                str(region.chromosome),
                region.start,
                region.end,
                strict_left=region.strict_left,
                strict_right=region.strict_right)

            if region.strand is not None:
                gene_strand = overlap['strand'].map({'+': 1, '-': -1})
                overlap = overlap.loc[gene_strand == region.strand]

            return overlap['gene_id']
        except KeyError:
            # TODO: print warning?
            return []

        return overlap


class Window(object):
    """Class respresenting a relative genomic window."""

    def __init__(self,
                 upstream,
                 downstream,
                 strand,
                 name=None,
                 strict_left=False,
                 strict_right=False):
        self.name = name
        self.upstream = upstream
        self.downstream = downstream
        self.strand = strand
        self.strict_left = strict_left
        self.strict_right = strict_right

    def apply(self, chromosome, position, strand):
        """Applies window to given position."""

        # Determine start/end position.
        if strand == 1:
            start = position - self.upstream
            end = position + self.downstream

            strict_left = self.strict_left
            strict_right = self.strict_right
        elif strand == -1:
            start = position - self.downstream
            end = position + self.upstream

            strict_right = self.strict_left
            strict_left = self.strict_right
        else:
            raise ValueError('Unknown value for strand ({})'.format(strand))

        strand = None if self.strand is None else self.strand * strand

        return Region(
            chromosome,
            start,
            end,
            strand,
            strict_left=strict_left,
            strict_right=strict_right)


Region = namedtuple('Region', [
    'chromosome', 'start', 'end', 'strand', 'strict_left', 'strict_right'
])
