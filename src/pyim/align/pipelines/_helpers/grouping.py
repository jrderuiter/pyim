import collections
import itertools
import operator

import heapq
import toolz

from collections import namedtuple


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

    def __str__(self):
        return 'PrioritySet(heap={}, set={})'\
            .format(str(self._heap), str(self._set))

    def __repr__(self):
        return str(self)


@toolz.curry
def groupby_reference(alignments, alignment_file=None):
    for reference, group in itertools.groupby(
            alignments, operator.attrgetter('reference_id')):
        if alignment_file is not None:
            reference = alignment_file.getrname(reference)
        yield reference, group


def groupby_position(alignments):
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
    aln_dict = collections.defaultdict(lambda: ([], []))

    curr_pos = 0
    for aln in alignments:
        # Check our ordering.
        if aln.reference_start < curr_pos:
            raise ValueError('Alignments not ordered by position')

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


GenomicPosition = namedtuple('GenomicPosition', 
                             ['chromosome', 'position', 'strand'])


def groupby_position_mate(alignments):
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
    aln_dict = collections.defaultdict(lambda: ([], []))

    # Only use proper pairs.
    alignments = (aln for aln in alignments if aln.is_proper_pair)

    # Limit ourselves to alignments from one chromosome (the first
    # encountered), as sort is only valid with the same chromosome.
    aln, alignments = toolz.peek(alignments)
    ref_name = aln.reference_name
    
    alignments = itertools.takewhile(
        lambda aln: aln.reference_name == ref_name, alignments)

    # We match position on the first pair. The second is stored until
    # needed and then returned together with the corresponding first pair.
    second_pairs = {}

    curr_pos = 0
    for aln in alignments:
        if aln.is_read2:
            second_pairs[aln.query_name] = aln
        else:
            # Check our ordering.
            if aln.reference_start < curr_pos:
                raise ValueError('Alignments not ordered by position')

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
                        fwd_mates = [second_pairs.pop(aln.query_name)
                                     for aln in fwd_grp]
                        fwd_pos = fwd_grp[0].reference_start
                        yield (GenomicPosition(ref_name, fwd_pos, 1), 
                               fwd_grp, fwd_mates)

                    if len(rev_grp) > 0:
                        rev_mates = [second_pairs.pop(aln.query_name)
                                     for aln in rev_grp]
                        rev_pos = rev_grp[0].reference_start
                        yield (GenomicPosition(ref_name, rev_pos, 1), 
                               rev_grp, rev_mates)
    
            except ValueError:
                pass

    # We're done, yield any remaining alignment groups.
    for _ in range(len(position_set)):
        fwd_grp, rev_grp = aln_dict.pop(position_set.pop())

        if len(fwd_grp) > 0:
            fwd_mates = [second_pairs.pop(aln.query_name) for aln in fwd_grp]
            fwd_pos = fwd_grp[0].reference_start
            yield (GenomicPosition(ref_name, fwd_pos, 1), fwd_grp, fwd_mates)

        if len(rev_grp) > 0:
            rev_mates = [second_pairs.pop(aln.query_name) for aln in rev_grp]
            rev_pos = rev_grp[0].reference_start
            yield (GenomicPosition(ref_name, rev_pos, 1), rev_grp, rev_mates)


@toolz.curry
def groupby_reference_position(alignments, alignment_file=None):
    chained = chain_groupby(
        alignments, [groupby_reference(alignment_file=alignment_file),
                     groupby_position])
    for res in chained:
        yield res


@toolz.curry
def groupby_barcode(alignments, barcode_map):
    # Group alignments by barcodes.
    groups = collections.defaultdict(list)
    for aln in alignments:
        barcode = barcode_map[aln.query_name]
        groups[barcode].append(aln)

    # Yield group together with barcode.
    for barcode, group in groups.items():
        yield barcode, group


def chain_groupby(iterable, groupby_funcs):
    grouped = groupby_funcs[0](iterable)

    if len(groupby_funcs) == 1:
        for key, group in grouped:
            if not isinstance(key, tuple):
                key = (key,)
            yield key, group
    else:
        for key, group in grouped:
            if not isinstance(key, tuple):
                key = (key,)
            for sub_key, sub_group in chain_groupby(group, groupby_funcs[1:]):
                yield key + sub_key, sub_group
