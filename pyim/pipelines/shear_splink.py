from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
from builtins import (ascii, bytes, chr, dict, filter, hex, input,
                      int, map, next, oct, open, pow, range, round,
                      str, super, zip)

from enum import Enum
from contextlib import contextmanager
from pathlib import Path

import numpy as np
import pandas as pd
import pysam
from skbio import io as skbio_io, DNASequence, SequenceCollection

from pyim.alignment.genome import Bowtie2Aligner
from pyim.alignment.vector import ExactAligner
from pyim.cluster import cluster_frame_merged

from .base import (Pipeline, GenomicExtractor,
                   InsertionIdentifier, genomic_distance)


class ShearSplinkPipeline(Pipeline):

    @classmethod
    def configure_argparser(cls, subparsers, name='shear_splink'):
        parser = subparsers.add_parser(name, help=name + ' help')

        parser.add_argument('input', type=Path)
        parser.add_argument('output', type=Path)
        parser.add_argument('reference', type=Path)
        parser.add_argument('transposon', type=Path)
        parser.add_argument('barcodes', type=Path)
        parser.add_argument('linker', type=Path)

        parser.add_argument('--contaminants', type=Path)
        parser.add_argument('--barcode_mapping', type=Path)
        parser.add_argument('--min_genomic_length', type=int, default=15)
        parser.add_argument('--threads', type=int, default=1)

        return parser

    @classmethod
    def from_args(cls, args):
        # Read transposon, barcode and linker sequences.
        transposon_seq = DNASequence.read(str(args['transposon']))

        linker_seq = DNASequence.read(str(args['linker']))

        barcode_seqs = SequenceCollection.read(
            str(args['barcodes']), constructor=DNASequence)

        # Read contaminants if supplied.
        contaminant_seqs = SequenceCollection.read(
            str(args['contaminants']), constructor=DNASequence) \
            if args['contaminants'] is not None else None

        # Read barcode map if supplied.
        if barcode_seqs is not None and args['barcode_mapping'] is not None:
            barcode_map = pd.read_csv(args['barcode_mapping'], sep='\t')
            barcode_map = dict(zip(barcode_map['barcode'],
                                   barcode_map['sample']))
        else:
            barcode_map = None

        # Setup extractor and identifier for pipeline.
        extractor = ShearSplinkExtractor(
            transposon_sequence=transposon_seq,
            transposon_aligner=ExactAligner(try_reverse=True),
            barcode_sequences=barcode_seqs,
            barcode_map=barcode_map,
            barcode_aligner=ExactAligner(try_reverse=False),
            linker_sequence=linker_seq,
            linker_aligner=ExactAligner(try_reverse=False),
            contaminant_sequences=contaminant_seqs,
            min_length=args['min_genomic_length'],
            threads=args['threads'])

        aligner = Bowtie2Aligner(args['reference'], bam_output=True)
        identifier = ShearSplinkIdentifier(min_mapq=0)

        return cls(extractor=extractor, aligner=aligner, identifier=identifier)


class ShearSplinkStatus(Enum):
    contaminant = 1
    no_transposon = 2
    no_linker = 3
    no_barcode = 4
    too_short = 5
    success = 6


class ShearSplinkExtractor(GenomicExtractor):

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
                barcode_aln = self._barcode_aligner.\
                    align_multiple(self._barcodes, read)

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
                        return ((genomic, barcode_aln.query_id),
                                self.STATUS.success)

    @classmethod
    def _read_input(cls, file_path, format_=None):
        format_ = 'fasta' if format_ is None else format_
        for read in skbio_io.read(str(file_path), format=format_,
                                  constructor=DNASequence):
            yield read

    @classmethod
    @contextmanager
    def _open_out(cls, file_path, format_=None):
        with file_path.open('w') as file_:
            yield file_

    @classmethod
    def _write_out(cls, sequence, fh, format_=None):
        format_ = 'fasta' if format_ is None else format_
        skbio_io.write(sequence, into=fh, format=format_)


class ShearSplinkIdentifier(InsertionIdentifier):

    def __init__(self, min_mapq=37, merge_distance=10):
        super().__init__()

        self._min_mapq = min_mapq
        self._merge_distance = merge_distance

    def identify(self, alignment_path, sample_map=None):
        bam_file = pysam.AlignmentFile(str(alignment_path), 'rb')

        # Collect insertions from alignments.
        insertions = []
        for ref_id in bam_file.references:
            alignments = bam_file.fetch(reference=ref_id)
            alignments = (aln for aln in alignments
                          if aln.mapping_quality >= self._min_mapq)

            # Group alignments by genomic position.
            aln_groups = self._group_alignments_by_position(
                alignments, barcode_map=sample_map)

            for (pos, strand, bc), alns in aln_groups:
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

        return insertions.sort(['seqname', 'location'])

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