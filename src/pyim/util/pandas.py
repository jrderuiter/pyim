from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
from builtins import (ascii, bytes, chr, dict, filter, hex, input,
                      int, map, next, oct, open, pow, range, round,
                      str, super, zip)


def reorder_columns(frame, order, drop_extra=False):
    if drop_extra:
        return frame[order]
    else:
        extra_cols = [c for c in frame.columns if c not in set(order)]
        return frame[list(order) + extra_cols]
