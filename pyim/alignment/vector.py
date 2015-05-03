from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
from builtins import (ascii, bytes, chr, dict, filter, hex, input,
                      int, map, next, oct, open, pow, range, round,
                      str, super, zip)

from skbio.alignment import local_pairwise_align_ssw


class VectorAlignment(object):

    def __init__(self, query_id, query_start, query_end, query_len,
                 target_id, target_start, target_end, target_strand,
                 target_len, type, identity, coverage):
        self.query_id = query_id
        self.query_start = query_start
        self.query_end = query_end
        self.query_len = query_len
        self.target_id = target_id
        self.target_start = target_start
        self.target_end = target_end
        self.target_strand = target_strand
        self.target_len = target_len
        self.type = type
        self.identity = identity
        self.coverage = coverage

    @property
    def score(self):
        return self.identity * self.coverage

    def reverse(self, read):
        read_len = len(read)

        return self.__class__(
            query_id=self.query_id,
            query_start=self.query_start,
            query_end=self.query_end,
            query_len=self.query_len,
            target_id=self.target_id,
            target_start=read_len - self.target_end,
            target_end=read_len - self.target_start,
            target_len=self.target_len,
            target_strand=1 if self.target_strand == -1 else 1,
            type=self.type, identity=self.identity, coverage=self.coverage
        )


class VectorAligner(object):

    def __init__(self, **kwargs):
        pass

    def align(self, query, target):
        raise NotImplementedError()

    def align_multiple(self, queries, target, how='unique'):
        alignments = filter(bool, (self.align(q, target) for q in queries))
        alignments = list(alignments)

        num_alignments = len(alignments)
        if num_alignments == 0:
            return None
        elif num_alignments == 1:
            return alignments[0]
        else:
            if how == 'unique':
                raise ValueError('Multiple matching queries for target {}'
                                 .format(target.id))
            elif how == 'any':
                return alignments[0]
            else:
                raise ValueError('Unknown value for how ({})'.format(how))


class ExactAligner(VectorAligner):

    def __init__(self, try_reverse=False):
        super().__init__()
        self._try_reverse = try_reverse

    def align(self, query, target):
        alignment = self._align_exact(query, target, query_ori=1)

        # Try reverse complement if first alignment failed.
        if alignment is None and self._try_reverse:
            alignment = self._align_exact(
                query.reverse_complement(), target, query_ori=-1)

        return alignment

    @staticmethod
    def _align_exact(query, target, query_ori):
        # Note that this alignment returns the first occurrence it finds,
        # later occurrences will not be found and are not checked for.
        try:
            index = target.sequence.index(query.sequence)
        except ValueError:
            return None
        else:
            q_len = len(query)

            return VectorAlignment(
                query_id=query.id, query_start=0, query_end=q_len,
                query_len=q_len, target_id=target.id, target_start=index,
                target_end=index + q_len, target_strand=query_ori,
                target_len=len(target), type='exact',
                identity=1.0, coverage=1.0)


class SswAligner(VectorAligner):

    def __init__(self, try_reverse=False, filters=None):
        super().__init__()
        self._try_reverse = try_reverse
        self._filters = filters

    def align(self, query, target):
        fwd_alignment = self._align_ssw(query, target, query_ori=1)

        if self._try_reverse:
            rev_alignment = self._align_ssw(
                query.reverse_complement(), target, query_ori=-1)

            if fwd_alignment is None:
                # Default to reverse if no forward.
                alignment = rev_alignment
            elif rev_alignment is None:
                # Default to forward if no reverse.
                alignment = fwd_alignment
            else:
                # Otherwise choose the best of the two.
                if rev_alignment.score > fwd_alignment.score:
                    alignment = rev_alignment
                else:
                    alignment = fwd_alignment
        else:
            alignment = fwd_alignment

        return alignment

    def _align_ssw(self, query, target, query_ori):
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

        aln = VectorAlignment(
            query_id=query.id, query_start=q_start, query_end=q_end,
            query_len=len(query), target_id=target.id, target_start=t_start,
            target_end=t_end, target_strand=query_ori, target_len=len(target),
            type='ssw', identity=identity, coverage=coverage)

        # Check if alignment passes any filter.
        if self._filters is None:
            return aln
        else:
            for filter_ in self._filters:
                if filter_(aln):
                    return aln
            return None


class ChainedAligner(VectorAligner):

    def __init__(self, aligners):
        super().__init__()
        self._aligners = aligners

    def align(self, query, target):
        aln = None

        for aligner in self._aligners:
            aln = aligner.align(query, target)
            if aln is not None:
                break

        return aln


def filter_identity(aln, min_identity):
    return aln.identity >= min_identity


def filter_score(aln, min_score):
    return aln.score >= min_score


def filter_end_match(aln, min_coverage=0.5, min_identity=1.0):
    return aln.target_end == aln.target_len and \
        aln.coverage >= min_coverage and aln.identity >= min_identity
