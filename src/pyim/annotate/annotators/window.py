from collections import namedtuple
from itertools import chain
from pathlib import Path

import pandas as pd

from pyim.util.pandas import GenomicDataFrame

from .base import Annotator, AnnotatorCommand, CisAnnotator
from ..util import filter_blacklist, select_closest, annotate_insertion


class WindowAnnotator(Annotator):
    """Window annotator class."""

    def __init__(self, genes, windows, closest=False, blacklist=None):
        super().__init__()

        self._windows = windows
        self._genes = genes

        self._closest = closest
        self._blacklist = blacklist

    @classmethod
    def from_window_size(cls, genes, window_size, **kwargs):
        """Creates instance using given window size."""

        window = Window(
            upstream=window_size,
            downstream=window_size,
            strand=None,
            name=None,
            strict_left=False,
            strict_right=False)

        return cls(genes=genes, windows=[window], **kwargs)

    def annotate(self, insertions):
        yield from chain.from_iterable((self._annotate_insertion(ins)
                                        for ins in insertions))

    def _annotate_insertion(self, insertion):
        # Identify overlapping features.
        hits = []
        for window in self._windows:
            region = window.apply(insertion.chromosome, insertion.position,
                                  insertion.strand)

            overlap = self._get_genes(region)
            overlap = overlap.assign(window=window.name)

            hits.append(overlap)

        hits = pd.concat(hits, axis=0, ignore_index=True)

        # Filter for closest/blacklist.
        if self._closest:
            hits = select_closest(insertion, hits)

        if self._blacklist is not None:
            hits = filter_blacklist(hits, self._blacklist)

        # Annotate insertion with identified hits.
        yield from annotate_insertion(insertion, hits)

    def _get_genes(self, region):
        try:
            overlap = self._genes.search(
                region.chromosome,
                region.start,
                region.end,
                strict_left=region.strict_left,
                strict_right=region.strict_right)

            if region.strand is not None:
                overlap = overlap.loc[overlap['strand'] == region.strand]
        except KeyError:
            overlap = GenomicDataFrame(pd.DataFrame().reindex(
                columns=self._genes.columns))

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


Region = namedtuple(
    'Region',
    ['chromosome', 'start', 'end', 'strand', 'strict_left', 'strict_right'])


class WindowAnnotatorCommand(AnnotatorCommand):
    """WindowAnnotator command."""

    name = 'window'

    def configure(self, parser):
        super().configure(parser)

        # Required arguments.
        parser.add_argument('--gtf', required=True, type=Path)

        # Optional arguments.
        parser.add_argument('--window_size', default=20000, type=int)

        parser.add_argument('--closest', default=False, action='store_true')
        parser.add_argument('--blacklist', nargs='+', default=None)

    def run(self, args):
        # Read insertions and genes.
        insertions = self._read_insertions(args.insertions)
        genes = self._read_genes_from_gtf(args.gtf)

        # Setup annotator.
        if args.cis_sites is not None:
            cis_sites = list(self._read_cis_sites(args.cis_sites))

            sub_annotator = WindowAnnotator.from_window_size(
                genes=genes,
                window_size=args.window_size,
                closest=args.closest,
                blacklist=args.blacklist)

            annotator = CisAnnotator(
                annotator=sub_annotator, genes=genes, cis_sites=cis_sites)
        else:
            annotator = WindowAnnotator.from_window_size(
                genes=genes,
                window_size=args.window_size,
                closest=args.closest,
                blacklist=args.blacklist)

        # Annotate insertions and write output.
        annotated = annotator.annotate(insertions)
        self._write_output(args.output, insertions=annotated)
