from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
from builtins import (ascii, bytes, chr, dict, filter, hex, input,
                      int, map, next, oct, open, pow, range, round,
                      str, super, zip)

from collections import namedtuple


class VectorAlignment(object):

    def __init__(self, query_id, query_start, query_end,
                 target_id, target_start, target_end, target_strand,
                 type, identity, coverage):
        self.query_id = query_id
        self.query_start = query_start
        self.query_end = query_end
        self.target_id = target_id
        self.target_start = target_start
        self.target_end = target_end
        self.target_strand = target_strand
        self.type = type
        self.identity = identity
        self.coverage = coverage

    @property
    def score(self):
        return self.identity * self.coverage

    def reverse(self, read):
        read_len = len(read)

        return self.__class__(
            query_id=self.query_id,
            query_start=self.query_start,
            query_end=self.query_end,
            target_id=self.target_id,
            target_start=read_len - self.target_end,
            target_end=read_len - self.target_start,
            target_strand=1 if self.target_strand == -1 else 1,
            type=self.type, identity=self.identity, coverage=self.coverage
        )
