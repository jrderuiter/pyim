import collections
import heapq
import itertools
import operator
from enum import Enum

import pysam
import numpy as np
import pandas as pd
from toolz import curry, map, pipe, merge_with
from toolz.curried import filter

import skbio
from skbio import DNA
from skbio.alignment import local_pairwise_align_ssw


# --- Model --- #

Alignment = collections.namedtuple(
    'Alignment',
    ['query_id', 'query_start', 'query_end', 'query_len',
     'target_id', 'target_start', 'target_end', 'target_len',
     'strand', 'identity', 'coverage', 'score'])

ExtractResult = collections.namedtuple(
    'ExtractResult', ['genomic_sequence', 'barcode', 'status'])


def reverse_alignment(aln):
    target_len = aln.target_len

    return Alignment(
        query_id=aln.query_id, query_start=aln.query_start,
        query_end=aln.query_end, query_len=aln.query_len,
        target_id=aln.target_id, target_start=target_len - aln.target_end,
        target_end=target_len - aln.target_start, target_len=target_len,
        strand=aln.strand * -1, type=aln.type, identity=aln.identity,
        coverage=aln.coverage, score=aln.score)


Insertion = collections.namedtuple(
    'Insertion', ['id', 'seqname', 'location',
                  'strand', 'sample', 'metadata'])


# --- Alignment --- #

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
def align_multiple(target, queries, align_func, **kwargs):
    alns = (align_func(target, query, **kwargs) for query in queries)
    alns = list(filter(bool, alns))

    if len(alns) == 0:
        return None
    elif len(alns) == 1:
        return alns[0]
    else:
        raise ValueError('Multiple alignments')


# --- Filtering --- #

def filter_alignment(alignment, filters):
    for filter_ in filters:
        if not filter_(alignment):
            return False
    return True


# --- Extract pipeline --- #

def extract(read):
    raise NotImplementedError()


def print_stats(results):
    # Iterate over results, counting statuses.
    status_counts = collections.defaultdict(int)

    for result in results:
        status_counts[result.status.name] += 1
        yield result

    # We're done, so print frequencies!
    total = sum(status_counts.values())
    for status, count in status_counts.items():
        percentage = (count / total) * 100
        print('{}: {} ({}%)'.format(status, count, percentage))


@curry
def write_sequences(results, file_path, format, mode='w',
                    compression='auto', compresslevel=9):
    """ Test docstring """
    with skbio.io.util.open(file_path, mode=mode, compression=compression,
                            compresslevel=compresslevel) as file_:
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


def consume(iterator, n=None):
    "Advance the iterator n-steps ahead. If n is none, consume entirely."
    # Use functions that consume iterators at C speed.
    if n is None:
        # Feed the entire iterator into a zero-length deque
        collections.deque(iterator, maxlen=0)
    else:
        # Advance to the empty slice starting at position n
        next(itertools.islice(iterator, n, n), None)


# --- Identify pipeline --- #

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


# --- ShearSplink --- #

class ShearSplinkStatus(Enum):
    contaminant = 1
    no_transposon = 2
    no_linker = 3
    no_barcode = 4
    multiple_barcodes = 5
    too_short = 6
    proper_read = 7


@curry
def shearsplink_extract(
        reads, transposon_sequence, barcode_sequences, linker_sequence,
        contaminant_sequences=None, transposon_func=None,
        barcode_func=None, linker_func=None, barcode_map=None):

    # Specify defaults for not provided aligners.
    if transposon_func is None:
        transposon_func = align_with_reverse(align_func=align_exact)

    if barcode_func is None:
        barcode_func = align_multiple(align_func=align_exact)

    if linker_func is None:
        linker_func = align_exact

    # Setup contaminant aligner if sequences are provided.
    if contaminant_sequences is not None:
        contaminant_func = align_multiple(queries=contaminant_sequences,
                                          align_func=align_exact)
    else:
        contaminant_func = None

    # Prime aligners with their respective sequences.
    transposon_func = transposon_func(query=transposon_sequence)
    barcode_func = barcode_func(queries=barcode_sequences)
    linker_func = linker_func(query=linker_sequence)

    # Extract and return results.
    extract_func = curry(_shearsplink_extract,
                         transposon_func=transposon_func,
                         barcode_func=barcode_func,
                         linker_func=linker_func,
                         contaminant_func=contaminant_func)

    for result in map(extract_func, reads):
        yield result


def _shearsplink_extract(
        read, transposon_func, barcode_func,
        linker_func, contaminant_func=None):
    """ Extracts the genomic sequence and barcode from the passed
        read. Reads containing contaminants are dropped. Reads are
        expected to look as follows:

            [barcode][transposon][genomic-sequence][linker]

        Each of these sequences is recognized by their corresponding
        alignment function. The barcode alignment identifies the
        barcode (and thus the sample) of the read, whilst the transposon
        and linker alignments are used to delineate the genomic sequence.

        The function returns an ExactResult tuple that contains the
        genomic sequence, barcode and a status flag. If any errors
        occur during the extraction, the genomic sequence and barcode
        values are None and the status flag indicates the underlying reason.
    """

    # Drop read if it contains a contaminant.
    if contaminant_func is not None and len(contaminant_func(read)) > 0:
        return ExtractResult(None, None, ShearSplinkStatus.contaminant)

    # Identify location of the transposon.
    transposon_aln = transposon_func(read)
    if transposon_aln is None:
        return ExtractResult(None, None, ShearSplinkStatus.no_transposon)

    # If transposon is on the reverse strand, flip the read and the
    # alignment to bring everything into the same (fwd) orientation.
    if transposon_aln.strand == -1:
        read = read.reverse_complement()
        transposon_aln = reverse_alignment(transposon_aln)

    # Identify barcode of the read.
    try:
        barcode_aln = barcode_func(read)
        if barcode_aln is None:
            return ExtractResult(None, None, ShearSplinkStatus.no_barcode)
    except ValueError:
        return ExtractResult(None, None, ShearSplinkStatus.multiple_barcodes)

    barcode = barcode_aln.query_id

    # Identify location of linker.
    linker_aln = linker_func(read)
    if linker_aln is None:
        return ExtractResult(None, None, ShearSplinkStatus.no_linker)

    # Extract genomic sequence using previous alignments.
    genomic = read[transposon_aln.target_end:linker_aln.target_start]

    return ExtractResult(genomic, barcode, ShearSplinkStatus.proper_read)


def shearsplink_identify(alignments):
    pass


def insertion_from_group(info, group):
    ref, pos, strand, bc = info

    # Get positions of the non-transposon ends of the alignment.
    end_field = 'reference_end' if strand == 1 else 'reference_start'
    end_positions = map(operator.attrgetter(end_field), group)

    # Calulate overall depth and unique end depth.
    depth = len(group)
    depth_unique = len(set(end_positions))

    metadata = dict(depth=depth, depth_unique=depth_unique)

    return Insertion(id=None, seq_name=ref, location=pos, strand=strand,
                     sample=bc, metadata=metadata)


def group_insertions(insertions, distance):
    # for insertion in insertions:
    #   check if we have an insertion from this sample in our collection
    #   if so, add to collection

    # - When did we last see SAMPLE_X?
    # - Which sample have we not seen within distance?
    pass


def merge_insertions(insertions):
    # Summarize location as mean.
    location = np.average([ins.location for ins in insertions])

    # Merge metadata by summing depths.
    metadata = merge_with(sum, *[ins.metadata for ins in insertions])

    # Take first insertion as reference for other attributes.
    ref = insertions[0]

    return Insertion(id=None, seqname=ref.seqname, location=location,
                     strand=ref.strand, sample=ref.sample, metadata=metadata)


# --- Main --- #

# Extraction.

seq1 = DNA('CACTGGCCACGCGAAGGTGC')
seq2 = DNA('GACCACTGGCCACGCGAAGG').reverse_complement()
seq3 = DNA('CGTTGGTCACTCTACCCACA')

transposon = DNA('TTTG', metadata=dict(id='transposon'))
barcodes = [DNA('AAAT', metadata=dict(id='BC01')),
            DNA('AAAA', metadata=dict(id='BC02'))]
linker = DNA('CCCG', metadata=dict(id='linker'))

reads = [DNA(str(barcodes[0]) + str(transposon) +
             str(seq1) + str(linker), metadata=dict(id='read_1')),
         DNA(str(transposon) + str(seq1) + str(linker))]

genomic_path = '/Users/Julian/Scratch/pyim/functional/genomic.fasta.gz'
barcode_path = '/Users/Julian/Scratch/pyim/functional/barcodes.txt'

barcode_map = pipe(
    reads,
    shearsplink_extract(transposon_sequence=transposon,
                        barcode_sequences=barcodes,
                        linker_sequence=linker),
    print_stats,
    filter(lambda r: r.status == ShearSplinkStatus.proper_read),
    filter(lambda r: len(r.genomic_sequence) >= 15),
    write_sequences(file_path=genomic_path, format='fasta',
                    compression='gzip', compresslevel=9),
    build_barcode_map)

barcode_frame = pd.DataFrame.from_records(
    iter(barcode_map.items()), columns=['read_id', 'barcode'])
barcode_frame.to_csv(barcode_path, sep='\t', index=False)


# Grouping.

bam = pysam.AlignmentFile('/Volumes/Datastore/Scratch/'
                          'lam-pcr-sjors/out/alignment.bam')
alns = itertools.islice(bam.fetch(), 0, 1000)

it = chain_groupby(
    itertools.islice(bam.fetch(multiple_iterators=True), 0, 1000),
    [curry(groupby_reference, alignment_file=bam), groupby_position])

barcode_map = collections.defaultdict(lambda: 'BC01')
it2 = chain_groupby(
    itertools.islice(bam.fetch(multiple_iterators=True), 0, 1000),
    [groupby_reference(alignment_file=bam),
     groupby_position,
     groupby_barcode(barcode_map=barcode_map)])
