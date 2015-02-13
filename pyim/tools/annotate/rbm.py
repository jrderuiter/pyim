__author__ = 'Julian'

import pandas as pd

from pyim_common.io import gff

from .base import Annotator


# Window format: (us, ua, ds, da)
WINDOW_PRESETS = {
    'SB': (20000, 10000, 25000, 5000),
    'MULV': (20000, 120000, 40000, 5000),
    'MMTV': (20000, 120000, 40000, 5000)
}


class RbmAnnotator(Annotator):

    def __init__(self, gtf_file, windows=None, preset=None,
                 homogeneity=None, feature_type='gene', id_column='id'):
        super(Annotator, self).__init__()

        if windows is None and preset is None:
            raise ValueError('Either windows or preset must be specified')

        if preset is not None:
            windows = WINDOW_PRESETS[preset]

        self._features = gff.read(gtf_file)
        self._feature_type = feature_type
        self._windows = windows
        self._homogeneity = homogeneity
        self._id_col = id_column

    def annotate(self, frame, feature_type=None, id_column=None):
        # Retrieve defaults set in advance if arguments not supplied.
        feature_type = self._feature_type if feature_type is None else feature_type
        id_column = self._id_col if id_column is None else id_column

        # Run RBM on frame, setting id_column as index.
        frame = frame.set_index(id_column, drop=False)
        result = rbm(frame, self._features, self._windows,
                     homogeneity=self._homogeneity, type_=feature_type)

        # Rename id column to correspond to input before returning.
        return result.rename(columns={'id': id_column})

    @classmethod
    def register_parser(cls, subparsers, name):
        parser = super(RbmAnnotator, cls).register_parser(subparsers, name)

        # Extra positional arguments.
        parser.add_argument('gtf_file')

        # Mutex arguments? (Should code properly)
        group = parser.add_mutually_exclusive_group(required=True)
        group.add_argument('--preset', default=None, choices=WINDOW_PRESETS.keys())
        group.add_argument('--windows', default=None, nargs='+', type=int)

        # Optional arguments.
        parser.add_argument('--feature_type', default='gene', choices=['gene', 'transcript'])
        parser.add_argument('--homogeneity', default=None, type=float)
        parser.add_argument('--id_column', default='id')

        return parser


def rbm(frame, features, windows, type_='gene', homogeneity=None):
    return pd.concat([_rbm_row_homogeneity(row, features, windows,
                                           type_=type_, homogeneity=homogeneity)
                      for _, row in frame.iterrows()], ignore_index=True)


def _rbm_row_homogeneity(row, features, windows, type_, homogeneity=None):
    hits = _rbm_row(row, features, windows, type_)

    # Include reverse orientation if strand is uncertain. In this case,
    # we mark the orientation as 'mixed', because the true ori is unclear.
    if homogeneity is not None and row.strand_homogeneity < homogeneity:
        row_rev = row.copy()
        row_rev.strand = 1 if row_rev.strand == -1 else -1
        hits_rev = _rbm_row(row_rev, features, windows, type_)

        hits = pd.concat([hits, hits_rev], ignore_index=True)
        hits['orientation'] = 'mixed'

    return hits.drop_duplicates()


def _rbm_row(row, features, windows, type_):
    us, ua, ds, da = windows
    rbm_windows = [('upstream', 'sense', us),
                   ('upstream', 'antisense', ua),
                   ('downstream', 'sense', ds),
                   ('downstream', 'antisense', da)]

    # Find hits in each of the RBM windows.
    hits = [_search_features(row, features, window, dir_, ori, type_)
            for dir_, ori, window in rbm_windows]
    hits = pd.concat(hits, ignore_index=True).drop_duplicates()

    # Summarize hits.
    hits = pd.concat([hits, gff.expand_ensembl_attrs(hits['attribute'])], axis=1)

    hit_rows = [{'id': row.name,
                 'gene_id': hit.gene_id,
                 'type': _hit_type(row, hit),
                 'orientation': _hit_orientation(row, hit)}
                for _, hit in hits.iterrows()]

    hit_frame = pd.DataFrame.from_records(
        hit_rows, columns=['id', 'gene_id', 'type', 'orientation'])

    return hit_frame


def _search_features(cis, features, window_size, direction, orientation, feature_type):
    features = features.ix[features['feature'] == feature_type]
    window = _search_window(cis, window_size, direction, orientation)
    return gff.find_overlap_frame(features, *window)


def _search_window(cis, window_size, direction='upstream', orientation='sense'):
    direction = -1 if direction == 'upstream' else 1
    window_direction = direction * cis.strand

    if window_direction == 1:
        start, end = cis.location, cis.location + window_size
    else:
        start, end = cis.location - window_size, cis.location

    if orientation == 'sense':
        feat_strand = cis.strand
    else:
        feat_strand = 1 if cis.strand == -1 else -1

    return cis.seqname, start, end, feat_strand


def _hit_type(ins, feature, features=None):
    if feature.feature == 'gene':
        return _hit_type_gene(ins, feature)
    elif feature.feature == 'transcript':
        return _hit_type_transcript(ins, feature, features)
    else:
        raise ValueError('Unsupported feature type {}'.format(feature.feature))


def _hit_type_gene(ins, gene):
    if ins.location < gene.start:
        type_ = 'upstream' if ins.strand == 1 else 'downstream'
    elif ins.location > gene.end:
        type_ = 'downstream' if ins.strand == 1 else 'upstream'
    else:
        type_ = 'intragenic'
    return type_


def _hit_type_transcript(ins, feature, features=None):
    type_ = _hit_type_gene(ins, feature)

    if type_ == 'intragenic' and features is not None:
        raise NotImplementedError('Detailed lookup not supported yet')

    return type_


def _hit_orientation(ins, feature):
    return 'sense' if ins.strand == feature.strand else 'antisense'