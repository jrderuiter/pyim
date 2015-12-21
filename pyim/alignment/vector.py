import collections

from skbio.alignment import local_pairwise_align_ssw
from toolz import curry


Alignment = collections.namedtuple(
    'Alignment',
    ['query_id', 'query_start', 'query_end', 'query_len',
     'target_id', 'target_start', 'target_end', 'target_len',
     'strand', 'identity', 'coverage', 'score'])


def reverse_alignment(aln):
    """Reverses strand of alignment object."""
    target_len = aln.target_len

    return Alignment(
        query_id=aln.query_id, query_start=aln.query_start,
        query_end=aln.query_end, query_len=aln.query_len,
        target_id=aln.target_id, target_start=target_len - aln.target_end,
        target_end=target_len - aln.target_start, target_len=target_len,
        strand=aln.strand * -1, type=aln.type, identity=aln.identity,
        coverage=aln.coverage, score=aln.score)


@curry
def align_exact(target, query, query_strand=1):
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
            strand=query_strand, identity=1.0, coverage=1.0, score=100)


@curry
def align_ssw(target, query, query_strand=1):
        ssw_aln = local_pairwise_align_ssw(target.sequence, query.sequence)

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
        identity = ssw_aln[0].fraction_same(ssw_aln[1])

        aln = Alignment(
            query_id=query.id, query_start=q_start, query_end=q_end,
            query_len=len(query), target_id=target.id, target_start=t_start,
            target_end=t_end, target_len=len(target), strand=query_strand,
            identity=identity, coverage=coverage,
            score=int(identity * coverage * 100))

        return aln


@curry
def align_with_reverse(target, query, align_func, query_strand=1, **kwargs):
    aln_fwd = align_func(target, query, query_strand=query_strand, **kwargs)
    aln_rev = align_func(target, query.reverse_complement(),
                         query_strand=query_strand * -1, **kwargs)

    if aln_fwd is None:
        return aln_rev
    elif aln_rev is None:
        return aln_fwd
    else:
        return aln_rev if aln_rev.score > aln_fwd.score else aln_fwd


@curry
def align_multiple(target, queries, align_func, return_first=False, **kwargs):
    alns = (align_func(target, query, **kwargs) for query in queries)
    alns = list(filter(bool, alns))

    if len(alns) == 0:
        return None
    elif len(alns) == 1 or return_first:
        return alns[0]
    else:
        raise ValueError('Multiple alignments')


# --- Filtering --- #

def filter_alignment(alignment, filters):
    for filter_ in filters:
        if not filter_(alignment):
            return False
    return True