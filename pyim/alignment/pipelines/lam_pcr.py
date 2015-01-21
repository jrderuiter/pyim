__author__ = 'Julian'

from numpy import sum as np_sum
from pysam import AlignmentFile

from pyim_common.model.insertion import Insertion

from .base import Pipeline
from ..util import alignment_frame


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
        alignment_path = input_path

        sam_file = AlignmentFile(alignment_path, 'rb')

        alignments = (aln for aln in sam_file if aln.mapping_quality >= self.min_mapq)

        # Map alignments to insertions and merge close-by insertions.
        insertions = self.map_insertions_chunked(sam_file, alignments, self.map_insertions)
        insertions = self.merge_insertions(insertions, self.merge_dist,
                                           metadata_func=self._merge_metadata)

        # Filter insertions with insufficient depth.
        insertions = [ins for ins in insertions if ins.metadata['depth'] >= self.min_depth]

        # Name insertions with a unique id.
        for i, ins in enumerate(insertions):
            ins.name = 'INS_{}'.format(i+1)

        Insertion.to_file(insertions, output_path)

    @staticmethod
    def map_insertions(sam_file, alignments):
        frame = alignment_frame(sam_file, alignments)

        # Define our anchor field to group on.
        is_fwd = frame['strand'] == 1
        frame.ix[is_fwd, 'anchor'] = frame.ix[is_fwd, 'start']
        frame.ix[~is_fwd, 'anchor'] = frame.ix[~is_fwd, 'end']

        groups = frame.groupby(['seqname', 'anchor', 'strand'])

        insertions = []
        for i, (_, group) in enumerate(groups):
            metadata = {
                'depth': group.shape[0],
                'depth_unique': group['anchor'].nunique(),
                'avg_mapq': group['mapq'].mean()
            }

            row = group.iloc[0]
            ins = Insertion(None, row.seqname, row.anchor, row.strand, None, metadata=metadata)

            insertions.append(ins)

        return insertions

    @staticmethod
    def merge_insertions(insertions, max_dist, metadata_func=None):
        return Insertion.merge_by_location(insertions, max_dist=max_dist,
                                           metadata_func=metadata_func)

    @staticmethod
    def _merge_metadata(ins_frame):
        total_depth = ins_frame['depth'].sum()

        avg_mapq = (ins_frame['depth'] / total_depth) * ins_frame['avg_mapq']
        avg_mapq = np_sum(avg_mapq)

        return {
            'depth_unique': ins_frame['depth_unique'].sum(),
            'depth': total_depth,
            'avg_mapq': avg_mapq
        }