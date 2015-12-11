from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
from builtins import (ascii, bytes, chr, dict, filter, hex, input,
                      int, map, next, oct, open, pow, range, round,
                      str, super, zip)


class Annotator(object):

    def __init__(self):
        super().__init__()

    @classmethod
    def configure_argparser(cls, subparsers, name='name'):
        raise NotImplementedError()

    @classmethod
    def from_args(cls, args):
        return cls(**args)

    def annotate(self, frame):
        raise NotImplementedError()


def get_closest(frame, id_col='insertion_id', distance_col='distance'):
    def _is_closest(x):
        abs_dist = x[distance_col].abs()
        return x.ix[abs_dist == abs_dist.min()]

    return (frame.groupby(id_col)
            .apply(_is_closest)
            .reset_index(drop=True))
