from collections import defaultdict
from itertools import groupby, chain

import numpy as np
import pysam

from pyim.model import Insertion
from pyim.vendor.frozendict import frozendict

from .external.cutadapt import cutadapt
from .external.bowtie2 import bowtie2


def trim_transposon(input_paths,
                    output_paths,
                    sequences,
                    threads=1,
                    verbose=False,
                    logger=None):
    """Trims input sequences for transposon sequence(s)."""

    # Check inputs.
    _check_reads(input_paths, output_paths)
    _check_adapters(sequences)

    # Assemble options.
    cutadapt_opts = dict(sequences)
    cutadapt_opts['--discard-untrimmed'] = True

    if len(input_paths) == 2:
        cutadapt_opts['--pair-filter=both'] = True

    return cutadapt(
        input_paths=input_paths,
        output_paths=output_paths,
        options=cutadapt_opts,
        threads=threads,
        verbose=verbose,
        logger=logger)


def trim_contaminants(input_paths,
                      output_paths,
                      sequences,
                      min_length=None,
                      threads=1,
                      verbose=False,
                      logger=None):
    """Trims contaminants from input sequences."""

    # Check inputs.
    _check_reads(input_paths, output_paths)
    _check_adapters(sequences)

    # Assemble options.
    cutadapt_opts = dict(sequences)

    if min_length is not None:
        cutadapt_opts['--minimum-length'] = min_length

    return cutadapt(
        input_paths=input_paths,
        output_paths=output_paths,
        options=cutadapt_opts,
        threads=threads,
        verbose=verbose,
        logger=logger)


def _check_reads(input_paths, output_paths=None):
    if not isinstance(input_paths, tuple):
        raise ValueError(
            'Input reads should be provided as a tuple of 1 (for single-end) '
            'or 2 (paired-end) file paths')

    if output_paths is not None:
        if not isinstance(output_paths, tuple):
            raise ValueError(
                'Output reads should be provided as a tuple of 1 (for '
                'single-end) or 2 (paired-end) file paths')

        if len(input_paths) != len(output_paths):
            raise ValueError('Number of input paths and output paths '
                             'should be the same')


def _check_adapters(adapters):
    """Checks adapters for valid keys."""

    if len(adapters) == 0:
        raise ValueError('No adapters given')

    valid_types = {'-a', '-b', '-g', '-A', '-B', '-G'}
    invalid_types = set(adapters.keys()) - valid_types

    if invalid_types:
        raise ValueError('Invalid adapter types: {}'.format(invalid_types))


def align_reads(read_paths,
                index_path,
                output_path,
                threads=1,
                verbose=False,
                logger=None):
    """Aligns reads to reference genome."""

    _check_reads(read_paths)

    bowtie2(
        read_paths=read_paths,
        index_path=index_path,
        output_path=output_path,
        extra_options={'--threads': threads},
        verbose=verbose,
        logger=logger)


def extract_insertions(alignment_path,
                       position_func,
                       sample_func,
                       min_mapq=30,
                       paired=False,
                       merge_dist=None,
                       min_support=0):
    """Extract insertions from alignments."""

    bam_file = pysam.AlignmentFile(str(alignment_path))

    try:
        summary = AlignmentSummary.from_alignments(
            iter(bam_file),
            position_func=position_func,
            sample_func=sample_func,
            min_mapq=min_mapq,
            paired=paired)
    finally:
        bam_file.close()

    if merge_dist is not None:
        summary = summary.merge_within_distance(merge_dist)

    yield from summary.to_insertions(min_support=min_support)


class AlignmentSummary(object):
    """Alignment summary class."""

    def __init__(self, values):
        self._values = values

    @property
    def values(self):
        """Returns alignment summary map."""
        return self._values

    @classmethod
    def from_alignments(cls,
                        alignments,
                        position_func,
                        sample_func,
                        paired=False,
                        min_mapq=30,
                        primary=True):
        """Constructs alignment summary from the given alignments."""

        # Optionally filter alignments.
        if primary:
            alignments = (aln for aln in alignments if not aln.is_secondary)

        if min_mapq is not None:
            alignments = (aln for aln in alignments
                          if aln.mapping_quality >= min_mapq)

        # Generate position/sample summaries.
        iter_func = cls._iter_paired if paired else cls._iter_single
        summaries = iter_func(alignments, position_func, sample_func)

        # Track alignment positions per sample. Each entry tracks positions
        # for a specific sample. Note that this dict contains two layers:
        #
        #  {'sample': {(ref, pos, strand): [linker_pos]}}
        #
        # The outer dict tracks positions for the given sample. The inner dict
        # (sample dict) actually tracks transposon positions as
        # (chrom, pos, strand) in the keys and linker positions (pos) in
        # the values, allowing us to derive transposon locations from the keys
        # and the number of 'unique ligation points' using the linker positions.
        alignment_map = defaultdict(lambda: defaultdict(list))
        for transposon_pos, linker_pos, sample in summaries:
            alignment_map[sample][transposon_pos].append(linker_pos)

        return cls(dict(alignment_map))

    @staticmethod
    def _iter_paired(alignments, position_func, sample_func):
        for mate1, mate2 in _iter_mates(alignments):
            transposon_pos, linker_pos = position_func(mate1, mate2)
            sample = sample_func(mate1, mate2)
            yield transposon_pos, linker_pos, sample

    @staticmethod
    def _iter_single(alignments, position_func, sample_func):
        for aln in alignments:
            transposon_pos, linker_pos = position_func(aln)
            sample = sample_func(aln)
            yield transposon_pos, linker_pos, sample

    def merge_within_distance(self, max_dist):
        """Merges summary values that are within given max dist."""

        merged_values = {}

        for sample, values in self._values.items():
            grouped = self._group_within_distance(values, max_dist=max_dist)
            new_values = dict(
                self._merge_sample_keys(grp, values) for grp in grouped)
            merged_values[sample] = new_values

        return self.__class__(merged_values)

    @classmethod
    def _group_within_distance(cls, sample_keys, max_dist):
        """Groups summary position keys that are within max dist."""

        # First, we sort by position and group by reference/strand.
        sorted_keys = sorted(sample_keys, key=lambda t: (t[2], t[0], t[1]))
        grouped_ref = groupby(sorted_keys, lambda t: (t[2], t[0]))

        # Then we identify and yield subgroups that are within max distance.
        grouped_pos = (cls._group_by_position(grp, max_dist=max_dist)
                       for _, grp in grouped_ref)  # yapf: disable

        yield from chain.from_iterable(grouped_pos)

    @staticmethod
    def _group_by_position(sample_keys, max_dist):
        sample_keys = list(sample_keys)

        curr_group = [sample_keys[0]]
        prev_pos = curr_group[0][1]

        for key in sample_keys[1:]:
            if (key[1] - prev_pos) > max_dist:
                yield curr_group
                curr_group = [key]
            else:
                curr_group.append(key)
            prev_pos = key[1]

        yield curr_group

    @staticmethod
    def _merge_sample_keys(key_group, sample_values):
        """Merges given summary position keys into a single entry."""

        # Get ref/strand.
        ref, _, strand = key_group[0]

        # Calculate (weighted) average position.
        grp_pos, grp_size = zip(*((k[1], len(sample_values[k]))
                                  for k in key_group))
        position = int(round(np.average(grp_pos, weights=grp_size)))

        # Determine new key and new entries.
        merged_key = (ref, position, strand)

        grp_values = (sample_values[key] for key in key_group)
        merged_value = list(chain.from_iterable(grp_values))

        return merged_key, merged_value

    def to_insertions(self, id_fmt='{sample}.INS_{num}', min_support=0):
        """Converts alignment map to a list of insertions."""

        i = 1

        for sample, sample_values in self._values.items():
            for key, ends in sample_values.items():
                ref, pos, strand = key
                support = len(set(ends))

                if support >= min_support:
                    metadata = frozendict(
                        depth=len(ends), depth_unique=support)

                    yield Insertion(
                        id=id_fmt.format(sample=sample, num=i),
                        sample=sample,
                        chromosome=ref,
                        position=pos,
                        strand=strand,
                        support=metadata['depth_unique'],
                        metadata=metadata)

                    i += 1


def _iter_mates(alignments):
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
