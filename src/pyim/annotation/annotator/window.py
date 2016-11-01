# pylint: disable=W0622,W0614,W0401
from __future__ import absolute_import, division, print_function
from builtins import *
# pylint: enable=W0622,W0614,W0401

import collections
import itertools

import toolz

from pyim.annotation import register_annotator
from pyim.util.tabix import GtfFile

# from ..filtering import filter_blacklist, select_closest
from ..util import build_interval_trees, numeric_strand


class WindowAnnotator(object):

    def __init__(self, reference_gtf, windows):
        if not isinstance(reference_gtf, GtfFile):
            reference_gtf = GtfFile(reference_gtf)

        self._windows = windows
        self._gtf = reference_gtf

        self._trees = None

    @classmethod
    def from_args(cls, args):
        window_size = args.window_size // 2
        windows = [Window(-window_size, window_size, strand=None,
                          name=None, strict_left=False, strict_right=False)]
        return cls(reference_gtf=args.reference_gtf, windows=windows)

    @classmethod
    def setup_args(cls, parser):
        # Required arguments.
        parser.add_argument('--reference_gtf', required=True)

        # Optional arguments.
        # parser.add_argument('--closest', default=False, action='store_true')
        parser.add_argument('--window_size', default=20000, type=int)

    def annotate(self, insertions):
        if self._trees is None:
            self._trees = build_interval_trees(self._gtf)

        queries = itertools.product(insertions, self._windows)
        annotated = itertools.chain.from_iterable(
            (self._annotate(ins, window, self._trees)
            for ins, window in queries))

        return annotated 

    def _annotate(self, ins, window, interval_trees):
        # Identify overlapping features.
        applied_window = window.apply(ins.chromosome, ins.position, ins.strand)
        features = list(applied_window.get_overlap(interval_trees))

        if len(features) > 0:
            for feature in features:
                feat_metadata = {'gene_id': feature['gene_id'],
                                 'gene_name': feature['gene_name']}

                if window.name is not None:
                    feat_metadata['window'] = window.name

                new_metadata = toolz.merge(ins.metadata, feat_metadata)

                yield ins._replace(metadata=new_metadata)
        else:
            yield ins


register_annotator('window', WindowAnnotator)


_Window =collections.namedtuple(
    'Window', ['start', 'end', 'strand', 'name',
               'strict_left', 'strict_right'])


class Window(_Window):
    __slots__ = ()

    def apply(self, chromosome, position, strand):
        # Determine start/end position.
        if strand == 1:
            start = position + self.start
            end = position + self.end

            strict_left = self.strict_left
            strict_right = self.strict_right
        elif strand == -1:
            start = position - self.end
            end = position - self.start

            strict_right = self.strict_left
            strict_left = self.strict_right
        else:
            raise ValueError('Unknown value for strand ({})'
                             .format(strand))

        # Determine new strand.
        if self.strand is not None:
            new_strand = self.strand * strand
        else:
            new_strand = None

        return AppliedWindow(chromosome, start, end, new_strand,
                             self.name, strict_left, strict_right)


_AppliedWindow = collections.namedtuple(
    'AppliedWindow', ['chromosome', 'start', 'end', 'strand',
                      'name', 'strict_left', 'strict_right'])


class AppliedWindow(_AppliedWindow):
    __slots__ = ()

    def get_overlap(self, interval_trees):
        # Find overlapping features.
        try:
            tree = interval_trees[self.chromosome]
            overlap = tree[self.start:self.end]
        except KeyError:
            overlap = []

        # Extract features.
        features = (interval[2] for interval in overlap)

        # Filter inclusive/exclusive if needed.
        if self.strict_left:
            features = (f for f in features if f['start'] > self.start)

        if self.strict_right:
            features = (f for f in features if f['end'] < self.end)

        # Filter for strand if needed.
        if self.strand is not None:
            features = (f for f in features
                        if numeric_strand(f['strand']) == self.strand)

        return features
