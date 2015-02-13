__author__ = 'Julian'

import pandas as pd

from pyim_common.io import gff

from pyim.tools.annotate.rbm import _hit_type, _hit_orientation
from .base import PandasDfWrapper


class InsertionFrame(PandasDfWrapper):

    @classmethod
    def read(cls, file_path):
        frame = pd.read_csv(file_path, sep='\t')
        frame = frame.rename(columns={'name': 'insertion_id'})
        return cls(frame)

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
            ref_genes = reference_gff.ix[reference_gff['feature'] == 'gene']
            ref_genes = pd.concat([ref_genes,
                                   gff.expand_ensembl_attrs(ref_genes['attribute'])], axis=1)
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
