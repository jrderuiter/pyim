__author__ = 'Julian'

import numpy as np
import pandas as pd
from pysam import AlignmentFile

from pyim.model.insertion import InsertionFrame

from .base import Pipeline, group_alignments_by_position


class LamPcrPipeline(Pipeline):
    # TODO: perform alignment within pipeline.

    def __init__(self, min_depth, merge_dist, min_mapq):
        super(LamPcrPipeline, self).__init__()

        self.min_depth = min_depth
        self.min_mapq = min_mapq
        self.merge_dist = merge_dist

    @classmethod
    def configure_argparser(cls, parser):
        parser = super(LamPcrPipeline, cls).configure_argparser(parser)

        parser.add_argument('--min-depth', default=5, type=int)
        parser.add_argument('--min-mapq', default=30, type=int)
        parser.add_argument('--merge-dist', default=10, type=int)

        return parser

    def run(self, input_path, output_path):
        insertions = self._insertions_from_alignments(input_path, min_mapq=self.min_mapq)

        # Filter insertions with insufficient depth.
        insertions = insertions.ix[insertions['depth_unique'] > self.min_depth]

        # Sort insertions by position and tag with id.
        insertions = insertions.sort(['seqname', 'location'])
        insertions['insertion_id'] = ['INS_{}'.format(i) for i in range(insertions.shape[0])]

        # Write insertions to disk.
        insertions.to_file(output_path)

    @staticmethod
    def _insertions_from_alignments(bam_path, min_mapq=30, merge_dist=10, bc_map=None):
        bam_file = AlignmentFile(bam_path, 'rb')

        # Collect insertions from alignments.
        insertions = []
        for ref_id in bam_file.references:
            print(ref_id)
            alignments = bam_file.fetch(reference=ref_id)
            alignments = (aln for aln in alignments if aln.mapping_quality >= min_mapq)

            aln_groups = group_alignments_by_position(alignments, barcode_map=bc_map)
            for (pos, strand, bc), alns in aln_groups:
                # Determine depth as the number of reads at this position.
                depth = len(alns)

                # Determine depth_unique by looking at differences in the
                # other position (end for fwd strand, start for rev strand).
                other_pos = (a.reference_end for a in alns) if strand == 1 else \
                    (a.reference_start for a in alns)
                other_pos = np.fromiter(other_pos, np.int, depth)
                depth_unique = len(np.unique(other_pos))

                insertions.append({'insertion_id': None, 'seqname': ref_id,
                                   'location': pos, 'strand': strand, 'sample': bc,
                                   'depth': depth, 'depth_unique': depth_unique})

        insertions = InsertionFrame(pd.DataFrame(insertions))

        # Merge insertions in close proximity.
        if merge_dist > 0:
            insertions = insertions.merge_by_location(
                max_dist=merge_dist, map_extra=[('depth', np.sum), ('depth_unique', np.sum)])

        return insertions
