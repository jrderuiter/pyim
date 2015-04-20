from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
from builtins import (bytes, dict, int, list, object, range, str,
                      ascii, chr, hex, input, next, oct, open,
                      pow, round, super, filter, map, zip)

import heapq


class PrioritySet(object):

    def __init__(self):
        self._heap = []
        self._set = set()

    def push(self, item, priority):
        if item not in self._set:
            heapq.heappush(self._heap, (priority, item))
            self._set.add(item)

    def pop(self):
        priority, item = heapq.heappop(self._heap)
        self._set.remove(item)
        return item

    def first(self):
        _, item = min(self._heap)
        return item

    def __len__(self):
        return len(self._heap)
