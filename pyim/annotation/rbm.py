from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
from builtins import (ascii, bytes, chr, dict, hex, input,
                      int, map, next, oct, open, pow, range, round,
                      str, super, zip)  # filter

from pathlib import Path

import pandas as pd

from toolz import curry, pipe, merge_with, keymap
from toolz.curried import filter, valfilter, valmap

from pyim.util.tabix import GtfFile
from pyim.util.pandas import reorder_columns

from .base import Annotator, get_closest
from .window import Window, apply_window, fetch_features, annotate_features


# Window format: (us, ua, ds, da)
WINDOW_SIZE_PRESETS = {
    'SB': (20000, 10000, 25000, 5000),
    'MULV': (20000, 120000, 40000, 5000),
    'MMTV': (20000, 120000, 40000, 5000)
}


class RbmAnnotator(Annotator):

    def __init__(self, gtf, window_sizes=None, preset=None,
                 feature_type='gene', closest=False, id_column='insertion_id'):
        super().__init__()

        if window_sizes is None:
            if preset is None:
                raise ValueError('Either windows or preset must be given')
            window_sizes = WINDOW_SIZE_PRESETS[preset]

        self._gtf = GtfFile(gtf)
        self._feature_type = feature_type

        self._closest = closest
        self._id_column = id_column

        self._windows = {
            'is': Window(None, 0, 1, 1, True, True),
            'ia': Window(None, 0, 1, -1, True, True),
            'us': Window(None, -window_sizes[0], 0, 1, True, False),
            'ua': Window(None, -window_sizes[1], 0, -1, True, False),
            'ds': Window(None, 1, window_sizes[2], 1, False, True),
            'da': Window(None, 1, window_sizes[3], -1, False, True)
        }

    @classmethod
    def configure_argparser(cls, subparsers, name='rbm'):
        parser = subparsers.add_parser(name, help=name + ' help')

        parser.add_argument('input', type=Path)
        parser.add_argument('output', type=Path)
        parser.add_argument('gtf')

        parser.add_argument('--feature_type', default='gene',
                            choices={'gene', 'transcript'})
        parser.add_argument('--id_column', default='insertion_id')
        parser.add_argument('--closest', default=False, action='store_true')

        group = parser.add_mutually_exclusive_group(required=True)
        group.add_argument('--preset')
        group.add_argument('--window_sizes', nargs=4, type=int)

        return parser

    def annotate(self, frame, type_='gene'):
        results = [self._annotate_row(row, self._windows,
                                      self._gtf, self._feature_type)
                   for _, row in frame.iterrows()]

        results = pd.concat(filter(lambda x: x is not None, results),
                            ignore_index=True)

        return results if not self._closest \
            else get_closest(results, id_col=self._id_column)

    @staticmethod
    def _annotate_row(row, windows, gtf, feature_type='gene'):
        strand = row.strand if hasattr(row, 'strand') else None

        # Fetch features for orientation, or for the forward orientation.
        apply_func = curry(apply_window, row.chrom,
                           row.position, strand or 1)
        windows_fwd = valmap(apply_func, windows)

        features = fetch_features_windows(
            gtf, windows_fwd, feature_type=feature_type)

        if strand is None:
            # Try again with reverse window orientation.
            apply_func = curry(apply_window, row.chrom, row.position, -1)
            windows_rev = valmap(apply_func, windows)

            features_rev = fetch_features_windows(
                gtf, windows_rev, feature_type=feature_type)

            # Reflect sense/antisense to match fwd windows.
            features_rev = keymap(
                curry(str_translate, table=str.maketrans('sa', 'as')),
                features_rev)

            features = merge_with(
                lambda frames: pd.merge(frames[0], frames[1], how='inner')
                    if len(frames) == 2 else frames[0],
                features, features_rev)

        annotated = {mech: annotate_features(row, features, mechanism=mech)
                     for mech, features in features.items()
                     if len(features) > 0}

        if len(annotated) > 0:
            frame = pd.concat(annotated.values(), ignore_index=True)
            return reorder_columns(frame, order=row.index)
        else:
            return None


def fetch_features_windows(gtf, windows, feature_type):
    return pipe(windows,
                valmap(lambda w: fetch_features(gtf, w, feature_type)),
                valfilter(lambda x: x is not None))


def str_translate(s, table):
    return str.translate(s, table)
