from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
from builtins import (ascii, bytes, chr, dict, filter, hex, input,
                      int, map, next, oct, open, pow, range, round,
                      str, super, zip)

from pathlib import Path

import numpy as np
import pandas as pd
import pysam

from pyim.alignment.genome import Bowtie2Aligner
from pyim.cluster import cluster_frame_merged
from pyim.io import read_fasta

from .base import (Pipeline, BasicGenomicExtractor,
                   InsertionIdentifier, genomic_distance)


class LamPcrPipeline(Pipeline):

    @classmethod
    def configure_argparser(cls, subparsers, name='lam_pcr'):
        parser = subparsers.add_parser('lam_pcr', help='lam_pcr help')

        parser.add_argument('transposon_sequence', type=Path)

        parser.add_argument('--barcode_sequences', type=Path)
        parser.add_argument('--barcode_mapping', type=Path)
        parser.add_argument('--linker_sequence', type=Path)

        return parser

    @classmethod
    def from_args(cls, args):
        transposon_seq = next(read_fasta(args['transposon_sequence']))

        # Read barcode and linker sequences if supplied.
        linker_seq = next(read_fasta(args['linker_sequence'])) \
            if args['linker_sequence'] is not None else None
        barcode_seqs = list(read_fasta(args['barcode_sequences'])) \
            if args['barcode_sequences'] is not None else None

        # Read barcode map if supplied.
        if barcode_seqs is not None and args['barcode_map'] is not None:
            barcode_map = pd.read_csv(args['barcode_map'], sep='\t')
            barcode_map = dict(zip(barcode_map['barcode'],
                                   barcode_map['sample']))
        else:
            barcode_map = None

        # Setup extractor and identifier for pipeline.
        extractor = BasicGenomicExtractor(
            transposon_seq, barcode_seqs, barcode_map, linker_seq)
        identifier = LamPcrIdentifier()

        return cls(extractor=extractor,
                   aligner=Bowtie2Aligner,
                   identifier=identifier)


class LamPcrIdentifier(InsertionIdentifier):

    def __init__(self, min_mapq=37, merge_distance=10):
        super().__init__()

        self._min_mapq = min_mapq
        self._merge_distance = merge_distance

    def identify(self, alignment_path, sample_map=None):
        bam_file = pysam.AlignmentFile(alignment_path, 'rb')

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