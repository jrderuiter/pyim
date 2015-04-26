from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
from builtins import (ascii, bytes, chr, dict, filter, hex, input,
                      int, map, next, oct, open, pow, range, round,
                      str, super, zip)

from collections import namedtuple


VectorAlignment = namedtuple(
    'Alignment',
    ['query_id', 'query_start', 'query_end',
     'target_id', 'target_start', 'target_end', 'target_strand',
     'type', 'identity', 'coverage']
)