import collections
import itertools
import operator

import skbio
from toolz import curry

from pyim.util import PrioritySet


def print_stats(results):
    # Iterate over results, counting statuses.
    status_counts = collections.defaultdict(int)

    for result in results:
        status_counts[result.status.name] += 1
        yield result

    # We're done, so print frequencies!
    print('\nExtract statistics:')

    total = sum(status_counts.values())
    for status, count in status_counts.items():
        percentage = (count / total) * 100
        print('{:>18}: {:>8} ({:05.2f}%)'.format(status, count, percentage))


@curry
def write_genomic_sequences(results, file_path, format='fastq',
                            mode='w', **io_kwargs):
    """ Test docstring """
    with skbio.io.open(file_path, mode, **io_kwargs) as file_:
        for result in results:
            skbio.io.write(result.genomic_sequence, into=file_, format=format)
            yield result


@curry
def build_barcode_map(results, sample_map=None):
    if sample_map is None:
        return {result.genomic_sequence.metadata['id']:
                result.barcode
                for result in results}
    else:
        return {result.genomic_sequence.metadata['id']:
                sample_map[result.barcode]
                for result in results}


@curry
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


@curry
def groupby_reference_position(alignments, alignment_file=None):
    chained = chain_groupby(
        alignments, [groupby_reference(alignment_file=alignment_file),
                     groupby_position])
    for res in chained:
        yield res


@curry
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
            for sub_key, sub_group in chain_groupby(group, groupby_funcs[1:]):
                yield key + sub_key, sub_group
