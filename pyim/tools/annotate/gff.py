__author__ = 'Julian'

import pandas
from pyim_common.io import gff
from pyim_common.util.frame import reorder_columns

from pyim.tools.annotate.base import Annotator


class GffAnnotator(Annotator):

    # TODO: Generalize gtf file type?
    def __init__(self, gtf_file, window):
        super(GffAnnotator, self).__init__()

        self.features = gff.read(gtf_file)
        self.window = window

    @classmethod
    def configure_argparser(cls, parser):
        parser = super(GffAnnotator, cls).configure_argparser(parser)
        parser.add_argument('--gtf-file', required=True)
        parser.add_argument('--window', default=20000, type=int)
        return parser

    def annotate_by_gene(self, insertions):
        features = self.features.ix[self.features['feature'] == 'gene']
        features = gff.expand_attr(features)

        window = self.window

        annotations = []
        for ins in insertions:
            range_ = (ins.seqname, ins.location - window, ins.location + window)
            overlap = gff.find_overlap(features, *range_)

            annotations += [self._to_dict(ins, feat) for _, feat in overlap.iterrows()]

        frame = pandas.DataFrame(annotations)
        frame.drop_duplicates(inplace=True)
        frame = reorder_columns(frame, ['insertion_id', 'gene_id', 'mechanism'])

        return frame

    def annotate_by_transcript(self, insertions):
        raise NotImplementedError

    @classmethod
    def _to_dict(cls, ins, feat):
        match_type, distance = _match_type(ins, feat)

        dict_ = {
            'insertion_id': ins.name,
            'gene_id': feat['gene_id'],
            'mechanism': match_type,
            'distance': distance
        }

        if 'gene_name' in feat:
            dict_['gene_name'] = feat['gene_name']

        return dict_


def _match_type(ins, feat):
    if feat.start <= ins.location <= feat.end:
        match, distance = 'within', 0
    else:
        if ins.location < feat.start:
            match = 'upstream' if feat.strand == 1 else 'downstream'
            distance = feat.start - ins.location
        else:
            match = 'downstream' if feat.strand == 1 else 'upstream'
            distance = ins.location - feat.end

    return match, distance
