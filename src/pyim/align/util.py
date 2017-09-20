from collections import defaultdict
from itertools import groupby, chain

import numpy as np

from pyim.model import Insertion
from pyim.vendor.frozendict import frozendict


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
        #   {'sample': {(ref, pos, strand): [linker_pos]}}
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
        for mate1, mate2 in iter_mates(alignments):
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
        sample_keys = iter(sample_keys)

        curr_group = [next(sample_keys)]
        prev_pos = curr_group[0][1]

        for key in sample_keys:
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

        for sample, sample_values in self._values.items():
            for i, (key, ends) in enumerate(sample_values.items()):
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
