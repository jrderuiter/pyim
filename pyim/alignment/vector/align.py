from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
from builtins import (ascii, bytes, chr, dict, filter, hex, input,
                      int, map, next, oct, open, pow, range, round,
                      str, super, zip)

from skbio.alignment import local_pairwise_align_ssw

from .model import VectorAlignment


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


class SswAligner(VectorAligner):

    def __init__(self, try_reverse=False):
        super().__init__()
        self._try_reverse = try_reverse

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

    @staticmethod
    def _align_ssw(query, target, query_ori):
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

        return VectorAlignment(
            query_id=query.id, query_start=q_start, query_end=q_end,
            target_id=target.id, target_start=t_start,
            target_end=t_end, target_strand=query_ori,
            type='ssw', identity=identity, coverage=coverage)


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
