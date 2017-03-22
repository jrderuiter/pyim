from collections import defaultdict
import itertools
import operator

import numpy as np
import toolz

from pyim.model import Insertion
from pyim.util.frozendict import frozendict


def extract_insertions(alignments,
                       func,
                       paired=False,
                       group_func=None,
                       merge_dist=None,
                       min_mapq=None,
                       min_support=1,
                       logger=None):
    """Extracts insertions from given alignments."""

    if logger is not None:
        logger.info('Summarizing alignments')
        logger.info('  %-18s: %s', 'Minimum mapq', min_mapq)

    alignments = filter_alignments(alignments, min_mapq=min_mapq)

    if group_func is None:
        summ_func = summarize_mates if paired else summarize_alignments
        summary = summ_func(alignments, func=func)

        if logger is not None:
            _log_insertions(logger, min_support, merge_dist)

        insertions = convert_summary_to_insertions(
            summary, merge_dist=merge_dist, min_support=min_support)
    else:
        if paired:
            raise NotImplementedError(
                'Grouping is not yet supported for paired-end data')
        else:
            aln_summaries = summarize_alignments_by_group(
                alignments, func, group_func=group_func)

            if logger is not None:
                _log_insertions(logger, min_support, merge_dist)

            insertion_grps = (
                convert_summary_to_insertions(
                    aln_summ,
                    min_support=min_support,
                    merge_dist=merge_dist,
                    sample=barcode,
                    id_fmt=barcode + '.INS_{}')
                for barcode, aln_summ in aln_summaries.items()) # yapf: disable

        insertions = list(itertools.chain.from_iterable(insertion_grps))

    return insertions


def _log_insertions(logger, min_support, merge_dist):
    logger.info('Converting to insertions')
    logger.info('  %-18s: %d', 'Minimum support', min_support)
    logger.info('  %-18s: %s', 'Merge distance', merge_dist)


def filter_alignments(alignments, only_primary=True, min_mapq=None):
    """Filters alignments on mapping quality, etc."""

    if only_primary:
        alignments = (aln for aln in alignments if not aln.is_secondary)

    if min_mapq is not None:
        alignments = (aln for aln in alignments
                      if aln.mapping_quality >= min_mapq)

    yield from alignments


def iter_mates(alignments):
    """Iterates over mate pairs in alignments."""

    cache = {}
    for aln in alignments:
        if aln.is_proper_pair:
            # Try to get mate from cache.
            mate = cache.pop(aln.query_name, None)

            if mate is None:
                # If not found, cache this mate.
                cache[aln.query_name] = aln
            else:
                # Otherwise, yield with mate.
                if aln.is_read1:
                    yield aln, mate
                else:
                    yield mate, aln


def summarize_alignments(alignments, func):
    """Summarizes alignments into a dict of chromosomal positions.

    This function summarizes an iterable of alignments into a dict that
    tracks the unique ends (ligation points) of the alignments for
    different genomic positions. The genomic positions are encoded as a tuple
    of (chromosome, position, strand) and are used as keys, whilst the
    ligation points are tracked as a list of positions.

    This dict is an intermediate used by other functions to derive insertions.

    Parameters
    ----------
    alignments : iterable[pysam.AlignedSegment]
        Alignments to summarize. May be prefiltered (on mapping quality
        for example), as this function does not perform any filtering itself.
    func : Function
        Function that takes an alignment and returns a Tuple containing
        (a) the location of the breakpoint with the transposon and (b) the
        position of the breakpoint with the linker (or the end of the read
        if no linkeris used). The former should be returned as a tuple of
        (chromosome, position, strand), whereas the latter should only be
        a position.

    Returns
    -------
    dict[(str, int, int), list[int]]
        Returns a dictionary mapping genomic positions, encoded as a
        (chromosome, position, strand) tuple to ligation points.

    """

    alignment_map = defaultdict(list)

    for aln in alignments:
        transposon_position, linker_position = func(aln)
        alignment_map[transposon_position].append(linker_position)

    return dict(alignment_map)


def summarize_mates(alignments, func):
    """Summarizes mate pairs into a dict of chromosomal positions."""

    alignment_map = defaultdict(list)

    for mate1, mate2 in iter_mates(alignments):
        transposon_position, linker_position = func(mate1, mate2)
        alignment_map[transposon_position].append(linker_position)

    return dict(alignment_map)


def summarize_alignments_by_group(alignments, func, group_func):
    """Summarizes groups of alignments into a dict of chromosomal positions."""

    # Take subgroups of alignments into account. This allows us to make
    # arbitrary subgroups of alignment summaries, for example by grouping
    # reads by sample barcodes.

    alignment_map = defaultdict(lambda: defaultdict(list))

    for aln in alignments:
        transposon_position, linker_position = func(aln)
        group = group_func(aln)

        if group is not None:
            alignment_map[group][transposon_position].append(linker_position)

    return {k: dict(v) for k, v in alignment_map.items()}


def merge_summary_within_distance(aln_summary, max_distance=10):
    """Merges alignment map entries that are within max_dist of each other."""

    grouped_keys = _groupby_position(aln_summary.keys(), max_distance)

    merged = dict(
        _merge_entries(aln_summary, key_grp) for key_grp in grouped_keys)

    return merged


def _groupby_position(alignment_keys, max_distance=10):
    """Groups alignment keys that are in close proximity for merging."""

    # First we sort by position and group by reference/strand.
    sorted_keys = sorted(alignment_keys, key=lambda t: (t[2], t[0], t[1]))
    grouped_keys = itertools.groupby(sorted_keys, lambda t: (t[2], t[0]))

    # Then we group the (position sorted) groups that are close together.
    grouped_pos = itertools.chain.from_iterable(
        _groupby_position_gen(
            v, max_distance=max_distance) for _, v in grouped_keys)

    return grouped_pos


def _groupby_position_gen(key_group, max_distance):
    key_iter = iter(key_group)

    prev = next(key_iter)
    curr_group = [prev]

    for key in key_iter:
        if (key[1] - prev[1]) <= max_distance:
            # Continue group.
            curr_group.append(key)
        else:
            # Start new group.
            yield curr_group
            curr_group = [key]

    yield curr_group


def _merge_entries(alignment_map, keys):
    # Calculate (weighted) average position.
    grp_pos, grp_size = zip(*((k[1], len(alignment_map[k])) for k in keys))
    pos = int(round(np.average(grp_pos, weights=grp_size)))

    # Generate new key/value.
    ref = keys[0][0]
    strand = keys[0][2]

    new_key = (ref, pos, strand)
    new_values = list(
        itertools.chain.from_iterable(alignment_map[k] for k in keys))

    return new_key, new_values


def convert_summary_to_insertions(aln_summary,
                                  min_support=1,
                                  merge_dist=None,
                                  id_fmt='INS_{}',
                                  **kwargs):
    """Converts an alignment map to a list of Insertions."""

    # Optionally merge insertions within x distance.
    if merge_dist is not None:
        aln_summary = merge_summary_within_distance(
            aln_summary, max_distance=merge_dist)

    # Convert to insertions.
    insertions = (_to_insertion(ref, pos, strand, ends, id_=None, **kwargs)
                  for i, ((ref, pos, strand), ends)
                  in enumerate(aln_summary.items())) # yapf: disable

    # Filter for support.
    insertions = (ins for ins in insertions if ins.support >= min_support)

    # Add ids.
    insertions = (ins._replace(id=id_fmt.format(i + 1))
                  for i, ins in enumerate(insertions))

    return list(insertions)


def _to_insertion(ref, pos, strand, ends, id_=None, **kwargs):
    metadata = toolz.merge({'depth': len(ends),
                            'depth_unique': len(set(ends))}, kwargs)
    return Insertion(
        id=id_,
        chromosome=ref,
        position=pos,
        strand=strand,
        support=metadata['depth_unique'],
        metadata=frozendict(metadata))
