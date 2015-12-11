from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
from builtins import (ascii, bytes, chr, dict, filter, hex, input,
                      int, map, next, oct, open, pow, range, round,
                      str, super, zip)

from collections import namedtuple
from functools import lru_cache
from pathlib import Path

import pandas as pd

from tkgeno.io import GtfFile
from tkgeno.util.pandas import reorder_columns

from .base import Annotator, get_closest


Window = namedtuple('Window', ['seqname', 'start', 'end', 'strand',
                               'incl_left', 'incl_right'])


def apply_window(seqname, location, strand, window):
    start = location + (window.start * strand)
    end = location + (window.end * strand)

    if strand == -1:
        start, end = end, start
        incl_left, incl_right = window.incl_right, window.incl_left
    else:
        incl_left, incl_right = window.incl_left, window.incl_right

    new_strand = strand * window.strand if window.strand is not None else None

    return Window(seqname, start, end, new_strand, incl_left, incl_right)


class WindowAnnotator(Annotator):

    def __init__(self, gtf_path, window_size, feature_type='gene',
                 id_column='insertion_id', closest=False):
        super().__init__()

        self._gtf = GtfFile(gtf_path)
        self._window = Window(
            seqname=None, start=-1 * window_size, end=window_size,
            strand=None, incl_left=True, incl_right=True)

        self._feature_type = feature_type
        self._closest = closest
        self._id_column = id_column

    @classmethod
    def configure_argparser(cls, subparsers, name='window'):
        parser = subparsers.add_parser(name, help=name + ' help')

        parser.add_argument('input', type=Path)
        parser.add_argument('output', type=Path)
        parser.add_argument('gtf')

        parser.add_argument('--feature_type', default='gene')
        parser.add_argument('--window_size', default=20000, type=int)

        parser.add_argument('--id_column', default='insertion_id')
        parser.add_argument('--closest', default=False, action='store_true')

        return parser

    def annotate(self, frame, type_='gene'):
        results = [self._annotate_row(row, self._window,
                                      self._gtf, self._feature_type)
                   for _, row in frame.iterrows()]

        results = pd.concat(filter(lambda x: x is not None, results),
                            ignore_index=True)

        return results if self._closest is not None \
            else get_closest(frame, id_col=self._id_column)

    @staticmethod
    def _annotate_row(row, window, gtf, feature_type='gene'):
        # Apply window for row.
        window = apply_window(row.seqname, row.location,
                              row.seqname, window)

        # Fetch features for row.
        features = fetch_features(gtf, window, feature_type=feature_type)

        # Annotate row with features, if any were found.
        if len(features) > 0:
            frame = annotate_features(row, features)
            return reorder_columns(frame, order=row.index)

        return None


@lru_cache(maxsize=64)
def fetch_features(gtf, window, feature_type):
    return gtf.get_region(feature=feature_type, **window._asdict())


def annotate_features(row, features, **kwargs):
    data = dict(row)
    data.update(dict(
        gene_id=features.gene_id,
        distance=[feature_distance(s, e, row.location)
                  for s, e in zip(features.start, features.end)]))
    data.update(**kwargs)

    return reorder_columns(pd.DataFrame(data), order=row.index)


def feature_distance(start, end, location):
    if start <= location <= end:
        return 0
    elif location > end:
        return location - end
    else:
        return location - start
