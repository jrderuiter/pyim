from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
# noinspection PyUnresolvedReferences
from builtins import (ascii, bytes, chr, dict, filter, hex, input,
                      int, map, next, oct, open, pow, range, round,
                      str, super, zip)
from future.utils import native_str

import logging
import pkg_resources
from collections import defaultdict
from multiprocessing import Pool

import pysam
import pandas as pd
import numpy as np
from scipy.spatial.distance import pdist
from skbio import io as skbio_io

# noinspection PyUnresolvedReferences
from pyim.io import fastq
from pyim.util import PrioritySet

logging.basicConfig(
    format='%(asctime)-15s %(message)s',
    datefmt='[%Y-%m-%d %H:%M:%S]',
    level=logging.INFO)


# --- Pipelines --- #

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

    def run(self, input_path, output_dir):
        logger = logging.getLogger()

        version = pkg_resources.get_distribution('pyim').version
        logger.info('--- PyIM v{} ---'.format(version))

        logger.info('Starting {} pipeline'.format(
            self.__class__.__name__.replace('Pipeline', '')))

        # Create directories if needed.
        if not output_dir.exists():
            output_dir.mkdir()

        # if input_path.suffix not in {'.bam', '.sam'}:
        #     genomic_path = output_dir / ('genomic' +
        #                                  ''.join(input_path.suffixes))
        #     barcode_path = output_dir / 'genomic.barcodes.txt'
        #
        #     # Extract genomic reads from input.
        #     logger.info('Extracting genomic sequences from reads')
        #
        #     _, barcodes = self._extractor.extract_file(
        #         input_path=input_path, output_path=genomic_path)
        #
        #     # Log statistics.
        #     total_reads = sum(self._extractor.stats.values())
        #
        #     logger.info('- Processed {} reads'.format(total_reads))
        #     logger.info('- Read statistics')
        #     for status in self._extractor.STATUS:
        #         count = self._extractor.stats[status]
        #         logger.info('\t- {}: {} ({:3.2f}%)'
        #                     .format(status.name, count,
        #                             (count / total_reads) * 100))
        #
        #     # Write out barcodes as frame.
        #     barcode_frame = pd.DataFrame.from_records(
        #         iter(barcodes.items()), columns=['read_id', 'barcode'])
        #     barcode_frame.to_csv(
        #         str(barcode_path), sep=native_str('\t'), index=False)
        #
        #     # Align to reference genome.
        #     logger.info('Aligning genomic sequences to reference')
        #     logger.info('- Using {} aligner (v{})'.format(
        #         self._aligner.__class__.__name__.replace('Aligner', ''),
        #         self._aligner.get_version()))
        #
        #     aln_path = self._aligner.align_file(
        #         file=genomic_path, output_dir=output_dir)
        # else:
        #     aln_path, barcodes = input_path, None

        aln_path = output_dir / 'alignment.bam'

        barcode_map = pd.read_csv(
            str(output_dir / 'genomic.barcodes.txt'), sep='\t')
        barcodes = dict(zip(barcode_map['read_id'], barcode_map['barcode']))

        # Identify transposon insertions.
        logger.info('Identifying insertions from alignment')

        insertions = self._identifier.identify(aln_path, barcode_map=barcodes)
        insertions.to_csv(str(output_dir / 'insertions.txt'),
                          sep=native_str('\t'), index=False)

        logger.info('--- Done! ---')


# --- Extractors --- #

# noinspection PyShadowingBuiltins
class GenomicExtractor(object):

    DEFAULT_IN_FORMAT = 'fasta'
    DEFAULT_OUT_FORMAT = 'fasta'

    def __init__(self, min_length=1, **kwargs):
        super().__init__()
        self._min_length = min_length
        self._stats = None

        self.reset_stats()

    @property
    def stats(self):
        return self._stats

    def reset_stats(self):
        self._stats = defaultdict(int)

    def extract(self, reads):
        for read in reads:
            result, status = self.extract_read(read)
            self._stats[status] += 1
            if result is not None:
                yield result

    def extract_read(self, read):
        raise NotImplementedError()

    def extract_from_file(self, file_path, format=None):
        format = self.DEFAULT_IN_FORMAT if format is None else format

        reads = skbio_io.read(str(file_path), format=format)
        for genomic, barcode in self.extract(reads):
                yield genomic, barcode

    def extract_to_file(self, reads, file_path, format=None):
        format = self.DEFAULT_OUT_FORMAT if format is None else format

        barcodes = {}
        with open(str(file_path), 'w') as file_:
            for genomic, barcode in self.extract(reads):
                barcodes[genomic.id] = barcode
                skbio_io.write(obj=genomic, format=format, into=file_)

        return file_path, barcodes

    def extract_file(self, input_path, output_path,
                     format_in=None, format_out=None):
        format_in = self.DEFAULT_IN_FORMAT if format_in is None else format_in
        format_out = self.DEFAULT_OUT_FORMAT \
            if format_out is None else format_out

        reads = skbio_io.read(str(input_path), format=format_in)
        return self.extract_to_file(reads, output_path, format=format_out)


class ParallelGenomicExtractor(GenomicExtractor):

    def __init__(self, min_length=1, threads=1, chunk_size=1000, **kwargs):
        super().__init__(min_length=min_length)

        self._threads = threads
        self._chunk_size = chunk_size

    def extract(self, reads):
        if self._threads == 1:
            for result in super().extract(reads):
                yield result
        else:
            pool = Pool(self._threads)

            for result, status in pool.imap_unordered(
                    self.extract_read, reads, chunksize=self._chunk_size):
                self._stats[status] += 1
                if result is not None:
                    yield result

            pool.close()
            pool.join()

    def extract_read(self, read):
        raise NotImplementedError()


# --- Identifiers --- #

class InsertionIdentifier(object):

    def __init__(self, **kwargs):
        super().__init__()

    def identify(self, alignment):
        raise NotImplementedError()

    @classmethod
    def _group_by_position_bam(cls, bam_path, barcode_map=None, min_mapq=0):
        bam_file = pysam.AlignmentFile(str(bam_path), 'rb')

        # Collect insertions from alignments.
        for ref_id in bam_file.references:
            alignments = bam_file.fetch(reference=ref_id)
            alignments = (aln for aln in alignments
                          if aln.mapping_quality >= min_mapq)

            # Group alignments by genomic position.
            aln_groups = cls._group_by_position_barcode(
                alignments, barcode_map=barcode_map)

            for (pos, strand, bc), alns in aln_groups:
                yield (ref_id, pos, strand, bc), alns

    @staticmethod
    def _group_by_position(alignments):
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

    @classmethod
    def _group_by_position_barcode(cls, alignments, barcode_map=None):
        grouped = cls._group_by_position(alignments)

        if barcode_map is None:
            for tup, grp in grouped:
                yield tup + (np.nan, ), grp
        else:
            for tup, grp in grouped:
                for bc, bc_grp in cls._split_by_barcode(grp, barcode_map):
                    yield tup + (bc, ), bc_grp

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
