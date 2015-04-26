from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
from builtins import (ascii, bytes, chr, dict, filter, hex, input,
                      int, map, next, oct, open, pow, range, round,
                      str, super, zip)

from collections import defaultdict

import numpy as np
from scipy.spatial.distance import pdist
from tkgeno.io import FastqFile

from pyim.alignment.vector.align import align_exact, align_multiple
from pyim.util import PrioritySet


class Pipeline(object):

    def __init__(self, extractor, aligner, identifier):
        super().__init__()
        self._extractor = extractor
        self._aligner = aligner
        self._identifier = identifier

    @classmethod
    def configure_argparser(cls, parser):
        raise NotImplementedError()

    @classmethod
    def from_args(cls, args):
        raise NotImplementedError()

    def run(self, input_path, output_path, work_dir):
        if input_path.suffix not in {'.bam', '.sam'}:
            # Assume Fastq for now.
            with FastqFile.open(input_path) as file_:
                # Extract genomic reads.
                genomic_path = work_dir / 'genomic.fastq.gz'
                self._extractor.extract_to_file(iter(file_), genomic_path)

            # Align to reference genome.
            aln_path = self._aligner.align_file(
                genomic_path, output_dir=work_dir)
        else:
            aln_path = input_path

        # Identify transposon insertions.
        insertions = self._identifier.identify(aln_path)
        insertions.to_csv(output_path, sep='\t')


class GenomicExtractor(object):

    def __init__(self, **kwargs):
        super().__init__()

    def extract(self, reads):
        raise NotImplementedError()

    def extract_to_file(self, reads, file_path):
        with FastqFile.open(file_path, mode='w') as file_:
            barcodes = {}
            for genomic, barcode in self.extract(reads):
                barcodes[genomic.id] = barcode
                FastqFile.write(file_, genomic)


class BasicGenomicExtractor(GenomicExtractor):

    def __init__(self, transposon_sequence, barcode_sequences=None,
                 barcode_map=None, linker_sequence=None):
        super().__init__()

        self._transposon_sequence = transposon_sequence
        self._barcodes = barcode_sequences
        self._barcode_map = barcode_map
        self._linker_sequence = linker_sequence

    def extract(self, reads):
        extract_func = self._dispatch()
        for read in reads:
            tup = self.extract_read(read, func=extract_func)
            if tup is not None:
                yield tup

    def extract_read(self, read, func=None):
        if func is None:
            func = self._dispatch()

        # Check for transposon sequence.
        transposon_aln = align_exact(self._transposon_sequence, read)

        # If we have a transposon sequence, try dispatching.
        if transposon_aln is not None:
            return func(read, transposon_aln)
        else:
            return None

    def _dispatch(self):
        has_barcodes = self._barcodes is not None
        has_linker = self._linker_sequence is not None

        if has_barcodes:
            return self._with_barcode_with_linker if has_linker \
                else self._with_barcode_no_linker
        else:
            return self._no_barcode_with_linker if has_linker \
                else self._no_barcode_no_linker

    def _with_barcode_with_linker(self, read, transposon_aln):
        # TODO: account for strand in extraction?
        # Match barcode and linker.
        barcode_aln = align_multiple(self._barcodes, read, align_exact)
        linker_aln = align_exact(self._linker_sequence, read)

        # Extract genomic sequence.
        if barcode_aln is not None and linker_aln is not None:
            genomic = read[barcode_aln.target_end:linker_aln.target_start]
            return genomic, barcode_aln.query_id
        else:
            return None

    def _with_barcode_no_linker(self, read, transposon_aln):
        barcode_aln = align_multiple(self._barcodes, read, align_exact)

        if barcode_aln is not None:
            genomic = read[barcode_aln.target_end:]
            return genomic, barcode_aln.query_id
        else:
            return None

    def _no_barcode_with_linker(self, read, transposon_aln):
        raise NotImplementedError()

    # noinspection PyMethodMayBeStatic
    def _no_barcode_no_linker(self, read, transposon_aln):
        return read[transposon_aln.target_end:], None


class InsertionIdentifier(object):

    def __init__(self, **kwargs):
        super().__init__()

    def identify(self, alignment):
        raise NotImplementedError()

    @classmethod
    def _group_alignments_by_position(cls, alignments, barcode_map=None):
        grouped = cls._group_alignments_by_pos(alignments)

        if barcode_map is None:
            for tup, grp in grouped:
                yield tup + (np.nan, ), grp
        else:
            for tup, grp in grouped:
                for bc, bc_grp in cls._split_by_barcode(grp, barcode_map):
                    yield tup + (bc, ), bc_grp

    @staticmethod
    def _group_alignments_by_pos(alignments):
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
        aln_dict = defaultdict(lambda: ([], []))

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

    @staticmethod
    def _split_by_barcode(alignments, barcode_map):
        split_groups = defaultdict(list)
        for aln in alignments:
            barcode = barcode_map[aln.query_name]
            split_groups[barcode].append(aln)

        for k, v in split_groups.items():
            yield k, v


def genomic_distance(insertions):
    loc = insertions['location']
    loc_2d = np.vstack([loc, np.zeros_like(loc)]).T
    dist = pdist(loc_2d, lambda u, v: np.abs(u-v).sum())
    return dist
