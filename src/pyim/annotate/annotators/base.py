import abc
from itertools import groupby, chain
import operator
from pathlib import Path

import numpy as np
import pandas as pd

from pyim.main import Command
from pyim.model import Insertion, CisSite
from pyim.vendor.genopandas import GenomicDataFrame

from ..util import annotate_insertion


class Annotator(abc.ABC):
    """Base annotator class."""

    @abc.abstractmethod
    def annotate(self, insertions):
        """Annotates given insertions with predicted target genes."""


class AnnotatorCommand(Command):
    """Base annotator command."""

    def configure(self, parser):
        parser.add_argument('--insertions', type=Path, required=True)
        parser.add_argument('--output', type=Path, required=True)

    @staticmethod
    def _read_insertions(insertion_path):
        return Insertion.from_csv(insertion_path, sep='\t')

    @staticmethod
    def _read_cis_sites(cis_path):
        return CisSite.from_csv(cis_path, sep='\t')

    @staticmethod
    def _read_genes_from_gtf(gtf_path):
        genes = GenomicDataFrame.from_gtf(
            gtf_path, filter_=lambda rec: rec['feature'] == 'gene')
        genes['strand'] = genes['strand'].map({'+': 1, '-': -1})
        return genes

    @staticmethod
    def _write_output(output_path, insertions):
        Insertion.to_csv(output_path, insertions, sep='\t', index=False)


class CisAnnotator(Annotator):
    """CIS annotator class."""

    def __init__(self, annotator, genes, cis_sites):
        super().__init__()

        if cis_sites is not None:
            # Copy unstranded CISs to both strands.
            cis_sites = self._expand_unstranded_sites(cis_sites)

        self._annotator = annotator
        self._genes = genes.set_index('gene_id', drop=False)
        self._cis_sites = cis_sites

    @staticmethod
    def _expand_unstranded_sites(cis_sites):
        """Copies unstranded CISs to both strands."""
        for cis in cis_sites:
            if np.isnan(cis.strand) or cis.strand is None:
                yield cis._replace(strand=1)
                yield cis._replace(strand=-1)
            else:
                yield cis

    def annotate(self, insertions):
        # Annotate cis sites.
        annotated_sites = self._annotator.annotate(self._cis_sites)
        cis_gene_mapping = self._extract_gene_mapping(annotated_sites)

        # Annotate insertions.
        annotated = chain.from_iterable((self._annotate_insertion(
            ins, cis_gene_mapping) for ins in insertions))

        yield from annotated

    @staticmethod
    def _extract_gene_mapping(annotated_sites):
        """Extracts CIS --> gene mapping from annotated CIS sites."""

        tuples = ((site.id, site.metadata['gene_id'])
                  for site in annotated_sites if 'gene_id' in site.metadata)
        grouped = groupby(sorted(tuples), key=operator.itemgetter(0))

        return {id_: set(tup[1] for tup in tuples) for id_, tuples in grouped}

    def _annotate_insertion(self, insertion, cis_gene_mapping):
        """Annotates an insertion using genes from mapping."""

        # Lookup hits in mapping.
        gene_ids = cis_gene_mapping.get(insertion.metadata['cis_id'], None)

        if gene_ids is None:
            hits = pd.DataFrame().reindex(columns=self._genes.columns)
        else:
            hits = self._genes.loc[gene_ids]

        # Annotate insertion with identified hits.
        yield from annotate_insertion(insertion, hits)
