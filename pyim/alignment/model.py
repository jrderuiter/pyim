
from collections import namedtuple

Alignment = namedtuple('Alignment', ['query_name', 'query_start', 'query_end', 'query_seq',
                                     'target_name', 'target_start', 'target_end', 'target_seq',
                                     'score', 'identity', 'alignment', 'type'])
