import collections

from skbio.alignment import local_pairwise_align_ssw
from toolz import curry


class Alignment(object):

    __slots__ = ('query_id', 'query_start', 'query_end', 'query_len',
                 'target_id', 'target_start', 'target_end', 'target_len',
                 'strand', 'identity', 'coverage', 'score', 'type')

    def __init__(self, query_id, query_start, query_end, query_len,
                 target_id, target_start, target_end, target_len,
                 strand, identity, coverage, type):
        self.query_id = query_id
        self.query_start = query_start
        self.query_end = query_end
        self.query_len = query_len

        self.target_id = target_id
        self.target_start = target_start
        self.target_end = target_end
        self.target_len = target_len

        self.strand = strand
        self.identity = identity
        self.coverage = coverage
        self.type = type

        self.score = int(identity * coverage * 100)

    def reverse(self):
        return Alignment(query_id=self.query_id,
                         query_start=self.query_start,
                         query_end=self.query_end,
                         query_len=self.query_len,
                         target_id=self.target_id,
                         target_start=self.target_len - self.target_end,
                         target_end=self.target_len - self.target_start,
                         target_len=self.target_len,
                         strand=self.strand * -1,
                         identity=self.identity,
                         coverage=self.coverage,
                         type=self.type)


@curry
def align_exact(target, query, query_strand=1):
    """Aligns query to target using exact matching."""

    # Note that this alignment returns the first occurrence it finds,
    # later occurrences will not be found and are not checked for.
    try:
        index = str(target).index(str(query))
    except ValueError:
        return None
    else:
        q_len = len(query)

        return Alignment(
            query_id=query.metadata.get('id', None), query_start=0,
            query_end=q_len, query_len=q_len,
            target_id=target.metadata.get('id', None), target_start=index,
            target_end=index + q_len, target_len=len(target),
            strand=query_strand, identity=1.0, coverage=1.0, type='exact')


@curry
def align_ssw(target, query, query_strand=1):
    """Aligns query to target using ssw aligner."""

    # Perform actual alignment.
    ssw_aln = local_pairwise_align_ssw(str(target), str(query))

    # Extract positions.
    pos = ssw_aln.start_end_positions()
    q_start, q_end = pos[1]
    t_start, t_end = pos[0]

    # Offset ends by one, making them exclusive
    # to match python conventions.
    q_end += 1
    t_end += 1

    # Calculate basic metrics.
    coverage = (q_end - q_start) / float(len(query))
    identity = ssw_aln[0].match_frequency(ssw_aln[1], relative=True)

    aln = Alignment(
        query_id=query.metadata.get('id', None), query_start=q_start,
        query_end=q_end, query_len=len(query),
        target_id=target.metadata.get('id', None), target_start=t_start,
        target_end=t_end, target_len=len(target), strand=query_strand,
        identity=identity, coverage=coverage, type='ssw')

    return aln


@curry
def align_with_reverse(target, query, align_func, query_strand=1, **kwargs):
    """Aligns query in both orientations to target sequence."""

    aln_fwd = align_func(target, query, query_strand=query_strand, **kwargs)
    aln_rev = align_func(target, query.reverse_complement(),
                         query_strand=query_strand * -1, **kwargs)
    return _pick_best(list(filter(bool, [aln_fwd, aln_rev])))


@curry
def align_multiple(target, queries, align_func, raise_error=False, **kwargs):
    """Aligns multiple queries to target sequence."""

    alignments = (align_func(target, query, **kwargs)
                  for query in queries)
    alignments = list(filter(bool, alignments))

    if len(alignments) > 1 and raise_error:
        raise ValueError('Multiple alignments')

    return _pick_best(alignments)


def _pick_best(alignments):
    """Picks best alignment from list (based on score)."""

    if len(alignments) == 0:
        return None
    if len(alignments) == 1:
        return alignments[0]
    else:
        best = alignments[0]
        for aln in alignments:
            if aln.score > best.score:
                best = aln
        return best


@curry
def align_chained(target, query, align_funcs, **kwargs):
    """Chains multiple vector alignment functions."""

    for func in align_funcs:
        aln = func(target, query, **kwargs)
        if aln is not None:
            return aln
    return None


def compose(align_func, try_reverse=False,
            filter='and', filters=None, **kwargs):
    """Helper function to build an aligner."""

    if try_reverse:
        align_func = align_with_reverse(align_func=align_func)

    if filters is not None:
        if filter == 'and':
            align_func = filter_and(align_func=align_func, filters=filters)
        elif filter == 'or':
            align_func = filter_or(align_func=align_func, filters=filters)
        else:
            raise ValueError('Filter should be either "or" or "and" (not {})'
                             .format(filter))

    return align_func(**kwargs)


# --- Filtering --- #

@curry
def filter_and(target, query, align_func, filters, **kwargs):
    """Performs AND of filters on resulting alignments."""

    alignment = align_func(target, query, **kwargs)
    for filter_ in filters:
        if not filter_(alignment):
            return None
    return alignment


@curry
def filter_or(target, query, align_func, filters, **kwargs):
    """Performs OR of filters on resulting alignments."""

    return not filter_and(target, query, align_func, filters, **kwargs)


@curry
def filter_score(alignment, min_score):
    """Checks if alignment has minimum score."""

    return alignment.score >= min_score


@curry
def filter_coverage(alignment, min_coverage, min_identity):
    """Checks if alignment is at end of read."""

    return ((alignment.coverage >= min_coverage) and
            (alignment.identity >= min_identity))

@curry
def filter_end_match(alignment):
    """Checks if alignment is at end of read."""

    return alignment.target_end == alignment.target_len
