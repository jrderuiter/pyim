from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
from builtins import (ascii, bytes, chr, dict, filter, hex, input,
                      int, map, next, oct, open, pow, range, round,
                      str, super, zip)


def identity_filter(aln, min_identity):
    return aln.identity >= min_identity


def coverage_filter(aln, min_coverage):
    return aln.coverage >= min_coverage
