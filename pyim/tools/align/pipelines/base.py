__author__ = 'Julian'

import heapq
from collections import defaultdict


class Pipeline(object):

    REQUIRED_OPTIONS = []

    def __init__(self, **kwargs):
        super(Pipeline, self).__init__()

    @classmethod
    def configure_argparser(cls, parser):
        return parser

#    @staticmethod
#    def map_insertions_chunked(sam_file, alignments, map_func):
#        chunks = chunk_alignments(alignments)
#        chunk_ins = [map_func(sam_file, c) for c in chunks]
#        return flatten_list(chunk_ins)


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


def group_alignments_by_position(alignments, barcode_map=None):
    grouped = _group_alignments_by_pos(alignments)

    if barcode_map is None:
        for tup, grp in grouped:
            yield tup + (None, ), grp
    else:
        for tup, grp in grouped:
            for bc, bc_grp in _split_by_barcode(grp, barcode_map):
                yield tup + (bc, ), bc_grp


def _group_alignments_by_pos(alignments):
    """ Groups alignments by their positions, grouping forward strand
        alignments with the same start position and reverse strand
        alignments with the same end position. Assumes alignments
        are all on a single reference sequence.
    """
    # Setup our collections for tracking reads and positions.
    #
    # The priority set is used to track positions with alignment groups,
    # ensuring that no position is listed twice (the set part) and
    # always giving the lowest position first (the priority part).
    #
    # The alignment dict contains two lists for each position with at
    # least one alignment, one for forward reads and one for reverse.
    # Any alignments encountered as position x in orientation o are added
    # to the corresponding entry dict[x][o] in the list, in which
    # o is encoded as {0,1}, with 1 being for reverse strand alignments.
    position_set = PrioritySet()
    aln_dict = defaultdict(lambda: ([], []))

    for aln in alignments:
        curr_pos = aln.reference_start

        # Add current read to collections.
        is_reverse = aln.is_reverse
        ref_pos = aln.reference_end if is_reverse else curr_pos
        aln_dict[ref_pos][bool(is_reverse)].append(aln)
        position_set.push(ref_pos, ref_pos)

        # Return any alignment groups before our current position.
        try:
            while position_set.first() < curr_pos:
                first_pos = position_set.pop()
                fwd_grp, rev_grp = aln_dict.pop(first_pos)
                if len(fwd_grp) > 0:
                    yield (fwd_grp[0].reference_start, 1), fwd_grp
                if len(rev_grp) > 0:
                    yield (rev_grp[0].reference_end, -1), rev_grp
        except ValueError:
            pass

    # We're done, yield any remaining alignment groups.
    for _ in range(len(position_set)):
        fwd_grp, rev_grp = aln_dict.pop(position_set.pop())
        if len(fwd_grp) > 0:
            yield (fwd_grp[0].reference_start, 1), fwd_grp
        if len(rev_grp) > 0:
            yield (rev_grp[0].reference_end, -1), rev_grp


def _split_by_barcode(alignments, barcode_map):
    split_groups = defaultdict(list)
    for aln in alignments:
        barcode = barcode_map[aln.query_name]
        split_groups[barcode].append(aln)

    for k, v in split_groups.items():
        yield k, v
