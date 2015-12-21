from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
from builtins import (ascii, bytes, chr, dict, filter, hex, input,
                      int, map, next, oct, open, pow, range, round,
                      str, super, zip)
from future.utils import native_str

from enum import Enum
from functools import partial
from pathlib import Path

import numpy as np
import pandas as pd
from skbio import DNA, SequenceCollection

from pyim.alignment.genome import Bowtie2Aligner
from pyim.alignment._vector import (ExactAligner, SswAligner, ChainedAligner,
                                    filter_score, filter_end_match)
from pyim.cluster import cluster_frame_merged

from ._base import (Pipeline, ParallelGenomicExtractor,
                    InsertionIdentifier, genomic_distance)


class ShearSplinkPipeline(Pipeline):

    @classmethod
    def configure_argparser(cls, subparsers, name='shear_splink'):
        parser = subparsers.add_parser(name, help=name + ' help')

        parser.add_argument('input', type=Path)
        parser.add_argument('output_dir', type=Path)
        parser.add_argument('reference', type=Path)
        parser.add_argument('transposon', type=Path)
        parser.add_argument('barcodes', type=Path)
        parser.add_argument('linker', type=Path)

        parser.add_argument('--contaminants', type=Path)
        parser.add_argument('--barcode_mapping', type=Path)
        parser.add_argument('--min_genomic_length', type=int, default=15)
        parser.add_argument('--min_depth', type=int, default=2)
        parser.add_argument('--min_mapq', type=int, default=37)

        parser.add_argument('--threads', type=int, default=1)

        return parser

    @classmethod
    def from_args(cls, args):
        # Read transposon, barcode and linker sequences.
        transposon_seq = DNA.read(str(args['transposon']))

        linker_seq = DNA.read(str(args['linker']))

        barcode_seqs = SequenceCollection.read(
            str(args['barcodes']), constructor=DNA)

        # Read contaminants if supplied.
        contaminant_seqs = SequenceCollection.read(
            str(args['contaminants']), constructor=DNA) \
            if args['contaminants'] is not None else None

        # Read barcode map if supplied.
        if barcode_seqs is not None and args['barcode_mapping'] is not None:
            barcode_map = pd.read_csv(str(args['barcode_mapping']),
                                      sep=native_str('\t'))
            barcode_map = dict(zip(barcode_map['barcode'],
                                   barcode_map['sample']))
        else:
            barcode_map = None

        # Setup transposon aligner.

        transposon_filters = [
            # Require at least 90% of the sequence to be matched.
            partial(filter_score, min_score=0.9)
        ]

        transposon_aligner = ChainedAligner(
            [ExactAligner(try_reverse=True),
             SswAligner(try_reverse=True, filters=transposon_filters)])

        # Setup linker aligner.
        linker_filters = [
            # Require at least 90% of the sequence to be matched.
            partial(filter_score, min_score=0.9),

            # Perfect match at the end of the read?
            partial(filter_end_match, min_coverage=0.5, min_identity=0.9)
        ]

        linker_aligner = ChainedAligner(
            [ExactAligner(try_reverse=False),
             SswAligner(try_reverse=False, filters=linker_filters)]
        )

        # Setup extractor and identifier for pipeline.
        extractor = ShearSplinkExtractor(
            transposon_sequence=transposon_seq,
            transposon_aligner=transposon_aligner,
            barcode_sequences=barcode_seqs,
            barcode_map=barcode_map,
            barcode_aligner=ExactAligner(try_reverse=False),
            linker_sequence=linker_seq,
            linker_aligner=linker_aligner,
            contaminant_sequences=contaminant_seqs,
            min_length=args['min_genomic_length'],
            threads=args['threads'])

        aligner = Bowtie2Aligner(args['reference'], bam_output=True,
                                 threads=args['threads'])
        identifier = ShearSplinkIdentifier(
            min_mapq=args['min_mapq'], min_depth=args['min_depth'])

        return cls(extractor=extractor, aligner=aligner, identifier=identifier)


class ShearSplinkStatus(Enum):
    contaminant = 1
    no_transposon = 2
    no_linker = 3
    no_barcode = 4
    multiple_barcodes = 5
    too_short = 6
    proper_read = 7


class ShearSplinkExtractor(ParallelGenomicExtractor):

    STATUS = ShearSplinkStatus

    def __init__(self, transposon_sequence, barcode_sequences, linker_sequence,
                 contaminant_sequences=None, transposon_aligner=None,
                 barcode_aligner=None, linker_aligner=None, barcode_map=None,
                 min_length=1, threads=1, chunk_size=1000):
        super().__init__(min_length=min_length,
                         threads=threads,
                         chunk_size=chunk_size)

        # Sequences.
        self._transposon = transposon_sequence
        self._contaminants = contaminant_sequences
        self._barcodes = barcode_sequences
        self._linker = linker_sequence

        # Aligners.
        self._transposon_aligner = transposon_aligner \
            if transposon_aligner is not None \
            else ExactAligner(try_reverse=True)

        self._contaminant_aligner = ExactAligner(try_reverse=True)

        self._barcode_aligner = barcode_aligner \
            if barcode_aligner is not None else ExactAligner()

        self._linker_aligner = linker_aligner \
            if linker_aligner is not None else ExactAligner()

        # Barcode map if given (maps barcodes to samples).
        self._barcode_map = barcode_map

    def extract_read(self, read):
        # Check for contaminants.
        if self._contaminants is not None:
            contaminant_aln = self._contaminant_aligner.\
                align_multiple(self._contaminants, read, how='any')

            if contaminant_aln is not None:
                return None, self.STATUS.contaminant

        # Check for transposon sequence.
        transposon_aln = self._transposon_aligner.align(
            self._transposon, read)

        if transposon_aln is None:
            # Missing transposon sequence.
            return None, self.STATUS.no_transposon
        else:
            # If we have a transposon sequence, continue.
            if transposon_aln.target_strand == -1:
                # If transposon is on the reverse strand, flip the
                # read and the alignment to bring everything downstream
                # into the same (fwd) orientation.
                read = read.reverse_complement()
                transposon_aln = transposon_aln.reverse(read)

            linker_aln = self._linker_aligner.align(
                self._linker, read)

            if linker_aln is None:
                # Missing linker sequence.
                return None, self.STATUS.no_linker
            else:
                try:
                    barcode_aln = self._barcode_aligner.\
                        align_multiple(self._barcodes, read)
                except ValueError:
                    return None, self.STATUS.multiple_barcodes

                if barcode_aln is None:
                    # Missing barcode sequence.
                    return None, self.STATUS.no_barcode
                else:
                    # Read is complete, return genomic part and barcode.
                    genomic = read[transposon_aln.target_end:
                                   linker_aln.target_start]

                    if len(genomic) < self._min_length:
                        return None, self.STATUS.too_short
                    else:
                        barcode = barcode_aln.query_id

                        if self._barcode_map is not None:
                            barcode = self._barcode_map[barcode]

                        return ((genomic, barcode),
                                self.STATUS.proper_read)


class ShearSplinkIdentifier(InsertionIdentifier):

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
        insertions = insertions.ix[
            insertions['depth_unique'] >= self._min_depth]

        # Add clonality annotation.
        insertions = self._annotate_clonality(insertions)

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
            # Check if merging is sane.
            assert len(set(frame['seqname'])) == 1
            assert len(set(frame['strand'])) == 1
            assert len(set(frame['sample'].astype(str))) == 1

            # Pick first row as reference for shared fields.
            ref = frame.iloc[0]

            # Calculate new location as mean, biased towards
            # insertions with more weight (a higher ULP).
            weighted_loc = np.average(frame.location,
                                      weights=frame['depth_unique'])
            weighted_loc = int(round(weighted_loc))

            return pd.Series(
                {'insertion_id': np.nan,
                 'seqname': ref['seqname'],
                 'location': weighted_loc,
                 'strand': ref['strand'],
                 'sample': ref['sample'],
                 'depth': frame['depth'].sum(),
                 'depth_unique': frame['depth_unique'].sum()},
                index=ref.index)

    @staticmethod
    def _annotate_clonality(ins_frame):
        groups = ins_frame.groupby('sample')

        if len(groups) > 0:
            clonality = groups.apply(lambda grp: grp['depth_unique'] /
                                     grp['depth_unique'].max())

            clonality.index = clonality.index.droplevel()
            clonality.name = 'clonality'
        else:
            clonality = pd.Series({'clonality': np.NaN},
                                  index=ins_frame.index)

        ins_frame_clonality = pd.concat([ins_frame, clonality], axis=1)

        return ins_frame_clonality
