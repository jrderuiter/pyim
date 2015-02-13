__author__ = 'Julian'

from collections import namedtuple


VectorAlignment = namedtuple(
    'Alignment',
    ['query_id', 'query_start', 'query_end',
     'target_id', 'target_start', 'target_end', 'target_strand',
     'type', 'identity', 'coverage']
)