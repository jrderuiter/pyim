__author__ = 'Julian'

from skbio.alignment import local_pairwise_align_ssw

from pyim.alignment.vector.model import VectorAlignment


def align_exact(query, target):
    alignment = _align_exact(query, target, query_ori=1)

    # Try reverse complement if first alignment failed.
    if alignment is None:
        alignment = _align_exact(
            query.reverse_complement(), target, query_ori=-1)

    return alignment


def _align_exact(query, target, query_ori=1):
    try:
        index = target.sequence.index(query.sequence)
    except ValueError:
        return None
    else:
        q_len = len(query)
        return VectorAlignment(
            query_id=query.id, query_start=0, query_end=q_len,
            target_id=target.id, target_start=index,
            target_end=index + q_len, target_strand=query_ori,
            type='exact', identity=1.0, coverage=1.0)


def align_ssw(query, target):
    alignment = _align_ssw(query, target, query_ori=1)

    # Try reverse complement if first alignment failed.
    if alignment is None:
        alignment = _align_ssw(query.reverse_complement(), target, query_ori=-1)

    return alignment


def _align_ssw(query, target, query_ori=1):
    ssw_aln = local_pairwise_align_ssw(target, query)

    # Extract positions.
    pos = ssw_aln.start_end_positions()
    q_start, q_end = pos[1]
    t_start, t_end = pos[0]

    # Offset end by one, making them exclusive to match python conventions.
    q_end += 1
    t_end += 1

    # Calculate basic metrics.
    coverage = (q_end - q_start) / float(len(query))
    identity = ssw_aln[0].fraction_same(ssw_aln[1])

    return VectorAlignment(
        query_id=query.id, query_start=q_start, query_end=q_end,
        target_id=target.id, target_start=t_start,
        target_end=t_end, target_strand=query_ori,
        type='ssw', identity=identity, coverage=coverage)


def align_chained(query, target, funcs):
    aln = None

    for func in funcs:
        aln = func(query, target)
        if aln is not None:
            break

    return aln


def align_multiple(queries, target, func,
                   resolve_func=None, warn=True, **kwargs):
    alignments = [func(q, target, **kwargs) for q in queries]
    alignments = [a for a in alignments if a is not None]

    alignment = None
    num_alignments = len(alignments)

    if num_alignments == 1:
        alignment = alignments[0]
    else:
        if num_alignments > 1:
            if resolve_func is not None:
                alignment = resolve_func(alignments)
            else:
                if warn:
                    print('WARNING: dropping multiple alignments for target {}'
                          .format(target.id))

    return alignment
