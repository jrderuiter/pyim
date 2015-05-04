from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
from builtins import (ascii, bytes, chr, dict, filter, hex, input,
                      int, map, next, oct, open, pow, range, round,
                      str, super, zip)

from contextlib import contextmanager
from enum import Enum
from pathlib import Path

import numpy as np
import pandas as pd
from skbio import DNASequence, SequenceCollection
from tkgeno.io import FastqFile

from pyim.alignment.genome import Bowtie2Aligner
from pyim.alignment.vector import ExactAligner
from pyim.cluster import cluster_frame_merged

from ._base import (Pipeline, FastqGenomicExtractor,
                    InsertionIdentifier, genomic_distance)


class LamPcrPipeline(Pipeline):

    @classmethod
    def configure_argparser(cls, subparsers, name='lampcr'):
        parser = subparsers.add_parser(name, help=name + ' help')

        parser.add_argument('input', type=Path)
        parser.add_argument('output_dir', type=Path)
        parser.add_argument('reference', type=Path)

        parser.add_argument('--contaminants', type=Path)
        parser.add_argument('--barcode_sequences', type=Path)
        parser.add_argument('--barcode_map', type=Path)
        parser.add_argument('--min_genomic_length', type=int, default=15)

        parser.add_argument('--min_depth', type=int, default=2)
        parser.add_argument('--min_mapq', type=int, default=37)

        parser.add_argument('--threads', type=int, default=1)

        return parser

    @classmethod
    def from_args(cls, args):

        # Read contaminant sequences.
        contaminant_seqs = SequenceCollection.read(
            str(args['contaminants']), constructor=DNASequence) \
            if args['contaminants'] is not None else None

        # Read barcode sequences if supplied.
        barcode_seqs = SequenceCollection.read(
            str(args['barcode_sequences']), constructor=DNASequence) \
            if args['barcode_sequences'] is not None else None

        # Read barcode map if supplied.
        if barcode_seqs is not None and args['barcode_map'] is not None:
            barcode_map = pd.read_csv(args['barcode_map'], sep='\t')
            barcode_map = dict(zip(barcode_map['barcode'],
                                   barcode_map['sample']))
        else:
            barcode_map = None

        # Setup extractor and identifier for pipeline.
        extractor = LamPcrExtractor(
            barcode_sequences=barcode_seqs,
            barcode_map=barcode_map,
            contaminant_sequences=contaminant_seqs,
            # threads=args['threads'],
            min_length=args['min_genomic_length'])

        aligner = Bowtie2Aligner(reference=args['reference'],
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
    too_short = 4
    proper_read = 5


class LamPcrExtractor(FastqGenomicExtractor):

    STATUS = LamPcrStatus

    def __init__(self, barcode_sequences=None, barcode_map=None,
                 contaminant_sequences=None, min_length=1,
                 threads=1, chunk_size=1000):
        super().__init__(min_length=min_length, threads=threads,
                         chunk_size=chunk_size)

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

        if self._barcodes is not None:
            try:
                bc_aln = self._barcode_aligner.align_multiple(
                    self._barcodes, read)
            except ValueError:
                return None, self.STATUS.duplicate_barcode

            if bc_aln is None:
                return None, self.STATUS.no_barcode
            else:
                # Extract genomic sequence.
                # TODO: check extraction for correctness.
                genomic = read[bc_aln.target_end:]

                # Lookup barcode.
                barcode = bc_aln.query_id
                if self._barcode_map is not None:
                    barcode = self._barcode_map[barcode]
        else:
            genomic, barcode = read, None

        if len(genomic) < self._min_length:
            return None, self.STATUS.too_short
        else:
            return (read, barcode), self.STATUS.proper_read

    @classmethod
    def _read_input(cls, file_path):
        with FastqFile.open(file_path) as file_:
            for i, read in enumerate(file_):
                if i % 100000 == 0 and i > 0:
                    print('Processed {} reads'.format(i))
                yield read

    @classmethod
    @contextmanager
    def _open_out(cls, file_path):
        if file_path.suffixes[-1] == '.gz':
            file_path = file_path.with_suffix('')

        with FastqFile.open(file_path, 'wt') as file_:
            yield file_


class LamPcrIdentifier(InsertionIdentifier):

    def __init__(self, min_depth=0, min_mapq=37, merge_distance=10):
        super().__init__()

        self._min_depth = min_depth
        self._min_mapq = min_mapq
        self._merge_distance = merge_distance

    def identify(self, alignment_path, barcode_map=None):
        insertions = []

        groups = self._group_by_position_bam(
            alignment_path, min_mapq=self._min_mapq)
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
