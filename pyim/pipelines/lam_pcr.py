from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
# noinspection PyUnresolvedReferences
from builtins import (ascii, bytes, chr, dict, filter, hex, input,
                      int, map, next, oct, open, pow, range, round,
                      str, super, zip)

from enum import Enum
from pathlib import Path

import numpy as np
import pandas as pd
from skbio import DNASequence, SequenceCollection

from pyim.alignment.genome import Bowtie2Aligner
from pyim.alignment.vector import ExactAligner
from pyim.cluster import cluster_frame_merged

from ._base import (Pipeline, GenomicExtractor,
                    InsertionIdentifier, genomic_distance)


class LamPcrPipeline(Pipeline):

    @classmethod
    def configure_argparser(cls, subparsers, name='lampcr'):
        parser = subparsers.add_parser(name, help=name + ' help')

        parser.add_argument('input', type=Path)
        parser.add_argument('output_dir', type=Path)
        parser.add_argument('reference', type=Path)

        parser.add_argument('--contaminants', type=Path, default=None)
        parser.add_argument('--transposon', type=Path, default=None)
        parser.add_argument('--barcodes', type=Path, default=None)
        parser.add_argument('--barcode_map', type=Path, default=None)
        parser.add_argument('--min_genomic_length', type=int, default=15)

        parser.add_argument('--min_depth', type=int, default=2)
        parser.add_argument('--min_mapq', type=int, default=37)

        parser.add_argument('--threads', type=int, default=1)

        return parser

    @classmethod
    def from_args(cls, args):

        # Read transposon sequence.
        transposon_seq = DNASequence.read(str(args['transposon'])) \
            if args['transposon'] is not None else None

        # Read contaminant sequences.
        contaminant_seqs = SequenceCollection.read(
            str(args['contaminants']), constructor=DNASequence) \
            if args['contaminants'] is not None else None

        # Read barcode sequences if supplied.
        barcode_seqs = SequenceCollection.read(
            str(args['barcodes']), constructor=DNASequence) \
            if args['barcodes'] is not None else None

        # Read barcode map if supplied.
        if barcode_seqs is not None and args['barcode_map'] is not None:
            barcode_map = pd.read_csv(str(args['barcode_map']), sep='\t')
            barcode_map = dict(zip(barcode_map['barcode'],
                                   barcode_map['sample']))
        else:
            barcode_map = None

        # Setup extractor and identifier for pipeline.
        extractor = LamPcrExtractor(
            transposon_sequence=transposon_seq,
            barcode_sequences=barcode_seqs,
            barcode_map=barcode_map,
            contaminant_sequences=contaminant_seqs,
            min_length=args['min_genomic_length'])

        aligner = Bowtie2Aligner(reference=args['reference'], bam_output=True,
                                 local=True, threads=args['threads'])

        identifier = LamPcrIdentifier(min_mapq=args['min_mapq'],
                                      min_depth=args['min_depth'])

        return cls(extractor=extractor,
                   aligner=aligner,
                   identifier=identifier)


class LamPcrStatus(Enum):
    contaminant = 1
    no_barcode = 2
    duplicate_barcode = 3
    no_transposon = 4
    too_short = 5
    proper_read = 6


class LamPcrExtractor(GenomicExtractor):

    DEFAULT_IN_FORMAT = 'fastq'
    DEFAULT_OUT_FORMAT = 'fastq'

    STATUS = LamPcrStatus

    def __init__(self, transposon_sequence=None,
                 barcode_sequences=None, barcode_map=None,
                 contaminant_sequences=None, min_length=1,
                 threads=1, chunk_size=1000):
        super().__init__(min_length=min_length, threads=threads,
                         chunk_size=chunk_size)

        self._transposon_sequence = transposon_sequence
        self._transposon_aligner = ExactAligner(try_reverse=False)

        self._barcode_aligner = ExactAligner(try_reverse=False)
        self._barcodes = barcode_sequences
        self._barcode_map = barcode_map

        self._contaminant_aligner = ExactAligner(try_reverse=True)
        self._contaminants = contaminant_sequences

    def extract_read(self, read):
        # Check for contaminants.
        if self._contaminants is not None:
            contaminant_aln = self._contaminant_aligner.\
                align_multiple(self._contaminants, read, how='any')

            if contaminant_aln is not None:
                return None, self.STATUS.contaminant

        # Check for a transposon sequence if specified.
        tr_aln = None

        if self._transposon_sequence is not None:
            tr_aln = self._transposon_aligner.align(
                self._transposon_sequence, read)

            if tr_aln is None:
                return None, self.STATUS.no_transposon

        # Check for barcode sequences if specified.
        bc_aln, barcode = None, None

        if self._barcodes is not None:
            try:
                bc_aln = self._barcode_aligner.align_multiple(
                    self._barcodes, read)
            except ValueError:
                return None, self.STATUS.duplicate_barcode

            if bc_aln is None:
                return None, self.STATUS.no_barcode

            # Lookup barcode.
            barcode = bc_aln.query_id
            if self._barcode_map is not None:
                barcode = self._barcode_map[barcode]

        # Extract the genomic sequence.
        if tr_aln is not None:
            genomic = read[tr_aln.target_end:]
        elif bc_aln is not None:
            genomic = read[bc_aln.target_end:]
        else:
            genomic = read

        # Check for minimum length.
        if len(genomic) < self._min_length:
            return None, self.STATUS.too_short

        # Return read, barcode and alignment status.
        return (read, barcode), self.STATUS.proper_read


class LamPcrIdentifier(InsertionIdentifier):

    def __init__(self, min_depth=0, min_mapq=37, merge_distance=10):
        super().__init__()

        self._min_depth = min_depth
        self._min_mapq = min_mapq
        self._merge_distance = merge_distance

    def identify(self, alignment_path, barcode_map=None):
        insertions = []

        groups = self._group_by_position_bam(
            alignment_path, min_mapq=self._min_mapq, barcode_map=barcode_map)
        for (ref_id, pos, strand, bc), alns in groups:
            # Determine depth as the number of reads at this position.
            depth = len(alns)

            # Determine depth_unique by looking at differences in the
            # other position (end for fwd strand, start for rev strand).
            other_pos = (a.reference_end for a in alns) if strand == 1 \
                else (a.reference_start for a in alns)
            depth_unique = len(set(other_pos))

            insertions.append(
                {'insertion_id': np.nan, 'seqname': ref_id,
                 'location': pos, 'strand': strand, 'sample': bc,
                 'depth': depth, 'depth_unique': depth_unique})

        # Create insertion frame.
        insertions = pd.DataFrame.from_records(
            insertions, columns=['insertion_id', 'seqname', 'location',
                                 'strand', 'sample', 'depth', 'depth_unique'])

        # Merge insertions in close proximity to account for sequencing errors.
        if self._merge_distance > 0:
            insertions = cluster_frame_merged(
                insertions, groupby=['seqname', 'sample', 'strand'],
                dist_func=genomic_distance, merge_func=self._merge_insertions,
                linkage='complete', criterion='distance',
                t=self._merge_distance)

        # Filter by min_depth.
        insertions = insertions.ix[insertions['depth_unique'] > self._min_depth]

        # Sort by coordinate and add identifiers.
        insertions = insertions.sort(['seqname', 'location'])

        insertions['insertion_id'] = ['INS_{}'.format(i+1)
                                      for i in range(len(insertions))]

        return insertions

    @classmethod
    def _merge_insertions(cls, frame):
        if len(frame) == 0:
            return frame.iloc[0]
        else:
            ref = frame.iloc[0]
            return pd.Series(
                {'insertion_id': np.nan,
                 'seqname': ref['seqname'],
                 'location': int(frame['location'].mean()),
                 'strand': ref['strand'],
                 'sample': ref['sample'],
                 'depth': ref['depth'].sum(),
                 'depth_unique': ref['depth_unique'].sum()},
                index=ref.index)
