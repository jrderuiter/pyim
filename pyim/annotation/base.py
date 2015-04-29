from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
from builtins import (ascii, bytes, chr, dict, filter, hex, input,
                      int, map, next, oct, open, pow, range, round,
                      str, super, zip)


class Annotator(object):

    def __init__(self, **kwargs):
        super().__init__()

    @classmethod
    def configure_argparser(cls, subparsers, name='name'):
        raise NotImplementedError()

    @classmethod
    def from_args(cls, args):
        raise NotImplementedError()

    def annotate(self, frame):
        raise NotImplementedError()
