__author__ = 'Julian'

from itertools import chain
from collections import OrderedDict

import numpy as np
import pandas as pd

from scipy.cluster.hierarchy import single, fcluster
from scipy.spatial.distance import pdist

from tkgeno.io.gff import GffFile

from pyim.tools.annotate.rbm import _hit_type, _hit_orientation
from .base import PandasDfWrapper


class InsertionFrame(PandasDfWrapper):

    REQUIRED_COLUMNS = ['insertion_id', 'seqname', 'location',
                        'strand', 'sample']

    def __init__(self, frame):
        # CHeck if all columns are present.
        for c in self.REQUIRED_COLUMNS:
            if c not in frame.columns:
                raise ValueError('Missing required column {}'.format(c))

        # Order columns.
        extra_cols = [c for c in frame.columns if c not in self.REQUIRED_COLUMNS]
        frame = frame[self.REQUIRED_COLUMNS + extra_cols]

        super(InsertionFrame, self).__init__(frame)

    @classmethod
    def read(cls, file_path):
        frame = pd.read_csv(file_path, sep='\t')
        frame = frame.rename(columns={'name': 'insertion_id'})
        return cls(frame)

    def to_file(self, file_path):
        self._frame.to_csv(file_path, sep='\t', index=False)

    def to_gff(self):
        gff_rows = (self._gff_row(row) for _, row in self._frame.iterrows())
        return GffFile(pd.DataFrame(gff_rows))

    @staticmethod
    def _gff_row(row, width=200):
        return {
            'seqname': row.seqname,
            'source': '.',
            'feature': 'insertion',
            'start': row.location - int(width / 2.0),
            'end': row.location + int(width / 2.0),
            'score': '.',
            'strand': row.strand,
            'frame': '.',
            'attribute': 'ID={id};NAME={id}'.format(id=row.insertion_id)
        }

    def annotate_with_clonality(self, depth_col, inplace=False):
        frame = self._frame if inplace else self._frame.copy()

        # Calculate clonality by dividing by highest depth.
        groups = frame.groupby('sample')
        frame['clonality'] = groups[depth_col].transform(lambda x: x / x.max())

        if not inplace:
            return self.__class__(frame)

    def annotate_by_cis_gene(self, cis_map, cis_genes, reference_gff=None, inplace=False):
        # First, we create a lookup from insertion --> gene.
        ins_genes = pd.merge(cis_map, cis_genes, on='cis_id')
        ins_genes = ins_genes[['insertion_id', 'gene_id']].drop_duplicates()

        # Now we join it to our insertion frame.
        merged = pd.merge(self._frame, ins_genes, on='insertion_id')

        # Finally, we annotate the gene hit type for each insertion.
        if reference_gff is not None:
            # Subset reference to genes.
            ref_genes = reference_gff.get_type('gene')
            ref_genes.parse_attribute('gene_id')
            ref_genes.set_index('gene_id', inplace=True)

            types = []
            for _, ins in merged.iterrows():
                feature = ref_genes.ix[ins.gene_id]
                types.append((_hit_type(ins, feature),
                              _hit_orientation(ins, feature)))

            merged['gene_type'], merged['gene_orientation'] = zip(*types)

        if inplace:
            self._frame = merged
        else:
            return self.__class__(merged)

    def merge_by_clustering(self, clusters, map_extra=None):
        merged = self._merge_by_clustering(self._frame, clusters, map_extra=map_extra)
        return self.__class__(merged)

    @classmethod
    def _merge_by_clustering(cls, frame, clusters, map_extra=None):
        # Group insertions by clusters.
        grouping = dict(zip(frame.index, clusters))
        grouped = frame.groupby(grouping)

        # Define our list of aggregation functions.
        map_funcs = [('insertion_id', _str_join),
                     ('seqname', 'first'),
                     ('location', np.mean),
                     ('strand', 'first'),
                     ('sample', 'first')]

        if map_extra is not None:
            map_funcs += map_extra

        # Actually aggregate the insertions.
        merged = grouped.agg(OrderedDict(map_funcs))

        return merged

    def merge_by_location(self, max_dist, map_extra=None, sort=True):
        grp_cols = ['seqname', 'strand']

        # Also group on sample if present.
        if any(self['sample'].notnull()):
            grp_cols += 'sample'

        # Merge insertions by sample, seqname, strand groups.
        # Groups are 'chunked' into insertions with local proximity
        # to avoid calculating a full NxN distance matrix.
        merged = []
        for _, grp in self._frame.groupby(by=grp_cols, as_index=False):
            for chunk in self._chunk_by_location(grp, max_dist=2 * max_dist):
                if len(chunk) > 1:
                    # Determine clusters using 1d distance measure.
                    dists = _dist_1d(chunk['location'])
                    clusters = fcluster(single(dists), t=max_dist, criterion='distance')

                    # Merge insertions using clustering.
                    merged_chunk = self._merge_by_clustering(
                        chunk, clusters, map_extra=map_extra)
                else:
                    merged_chunk = chunk

                merged.append(merged_chunk)

        # Concat into single dataframe.
        merged = pd.concat(merged, ignore_index=True)
        if sort:
            merged.sort(['seqname', 'strand'], inplace=True)

        return self.__class__(merged)

    @staticmethod
    def _chunk_by_location(group, max_dist):
        group = group.sort('location')

        chunks = np.cumsum(np.diff(group['location']) > max_dist)
        chunk_dict = dict(zip(group.index, chain([0], chunks)))

        for _, chunk in group.groupby(chunk_dict):
            yield chunk


def _str_join(values, sep=';'):
    values = [v for v in values if v is not None]
    if len(values) > 0:
        return sep.join(values)
    else:
        return None


def _dist_1d(loc):
    loc_2d = np.vstack([loc, np.zeros_like(loc)]).T
    dist = pdist(loc_2d, lambda u, v: np.abs(u-v).sum())
    return dist
