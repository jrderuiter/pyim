
import pandas as pd
from collections import namedtuple

Alignment = namedtuple('Alignment', ['query_name', 'query_start', 'query_end', 'query_seq',
                                     'target_name', 'target_start', 'target_end', 'target_seq',
                                     'score', 'identity', 'alignment', 'type'])

def alignments_to_frame(alns):
    if len(alns) == 0: return None
    values = alns.values() if type(alns) == dict else alns
    return pd.DataFrame(values, columns=values[0]._fields)

#class Alignment(object):
#
#    def __init__(self, query_name, query_start, query_end, query_seq,
#                 target_name, target_start, target_end, target_seq,
#                 score, identity, alignment, type):
#        self.query_name = query_name
#        self.query_start = query_start
#        self.query_end = query_end
#        self.query_seq = query_seq
#
#        self.target_name = target_name
#        self.target_start = target_start
#        self.target_end = target_end
#        self.target_seq = target_seq
#
#        self.score = score
#        self.identity = identity
#        self.alignment = alignment
#        self.type = type
#
#    def __len__(self):
#        return self.query_end - self.query_start
#
#    def aligned_seq_query(self):
#        return self.query_seq[self.query_start:self.query_end]
#
#    def aligned_seq_target(self):
#        return self.target_seq[self.target_start:self.target_end]
#