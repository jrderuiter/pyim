from collections import defaultdict
import itertools
import logging
import operator

from frozendict import frozendict
import numpy as np
import pysam
import toolz

from pyim.model import Insertion


def fetch_alignments(bam_path, only_primary=True, min_mapq=None):
    bam_file = pysam.AlignmentFile(str(bam_path))

    try:
        alignments = iter(bam_file)

        if only_primary:
            alignments = (aln for aln in alignments if not aln.is_secondary)

        if min_mapq is not None:
            alignments = (aln for aln in alignments
                          if aln.mapping_quality >= min_mapq)

        yield from alignments
    finally:
        bam_file.close()


def summarize_alignments(alignments):
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

    Returns
    -------
    dict[(str, int, int), list[int]]
        Returns a dictionary mapping genomic positions, encoded as a
        (chromosome, position, strand) tuple to ligation points.

    """
    alignment_map = defaultdict(list)

    for aln in alignments:
        tup = _process_alignment(aln)
        if tup is not None:
            alignment_map[tup[0]].append(tup[1])

    return dict(alignment_map)


def summarize_alignments_by_group(alignments, group_func):
    # Take subgroups of alignments into account. This allows us to make
    # arbitrary subgroups of alignment summaries, for example by grouping
    # reads by sample barcodes.
    alignment_map = defaultdict(lambda: defaultdict(list))

    for aln in alignments:
        tup = _process_alignment(aln)
        if tup is not None:
            grp = group_func(aln)
            if grp is not None:
                alignment_map[grp][tup[0]].append(tup[1])

    return {k: dict(v) for k, v in alignment_map.items()}


def _process_alignment(aln):
    if aln.reference_id != -1:
        ref = aln.reference_name

        if aln.is_reverse:
            transposon_pos = aln.reference_end
            linker_pos = aln.reference_start
            strand = -1
        else:
            transposon_pos = aln.reference_start
            linker_pos = aln.reference_end
            strand = 1

        key = (ref, transposon_pos, strand)

        return key, linker_pos
    else:
        return None


def extract_barcode_mapping(reads, barcodes, barcode_mapping=None):

    # Create barcode/sample dict.
    barcode_dict = {bc.name: bc.sequence for bc in barcodes}

    if barcode_mapping is not None:
        barcode_dict = {sample: barcode_dict[barcode]
                        for barcode, sample in barcode_mapping.items()}

    # Build mapping.
    mapping = {}

    for read in reads:
        # Check each barcode for match in read.
        matched = [k for k, v in barcode_dict.items() if v in read.sequence]

        if len(matched) == 1:
            # Record single matches.
            name = read.name.split()[0]
            mapping[name] = matched[0]
        elif len(matched) > 1:
            logging.warning('Skipping %s due to multiple matching barcodes',
                            read.name.split()[0])

    return mapping


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
                                  merge_distance=None,
                                  id_fmt='INS_{}',
                                  **kwargs):
    """Converts an alignment map to a list of Insertions."""

    # Optionally merge insertions within x distance.
    if merge_distance is not None:
        aln_summary = merge_summary_within_distance(
            aln_summary, max_distance=merge_distance)

    # Convert to insertions.
    insertions = (_to_insertion(ref, pos, strand, ends, id_=None, **kwargs)
                  for i, ((ref, pos, strand), ends)
                  in enumerate(aln_summary.items())) # yapf: disable

    # Filter for support.
    insertions = (ins for ins in insertions if ins.support >= min_support)

    # Sort by depth and add IDs.
    insertions = sorted(insertions, key=operator.attrgetter('support'))[::-1]
    insertions = [ins._replace(id=id_fmt.format(i + 1))
                  for i, ins in enumerate(insertions)]

    return insertions


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
