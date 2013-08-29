

import pandas
from collections import namedtuple


Alignment = namedtuple('Alignment', ['query_name', 'query_start', 'query_end', 'query_seq',
                                     'target_name', 'target_start', 'target_end', 'target_seq',
                                     'score', 'identity', 'alignment', 'type'])
def alignments_to_frame(alns):
    if len(alns) == 0: return None
    values = alns.values() if type(alns) == dict else alns
    return pandas.DataFrame(values, columns=values[0]._fields)
