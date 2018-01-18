import abc

import pandas as pd

from pyim.model import InsertionSet  #, CisSite

from ..util import annotate_matches, filter_matches


class Annotator(abc.ABC):
    """Base annotator class."""

    @abc.abstractmethod
    def match(self, insertions, genes):
        """"Matches insertion to predicted target genes."""

    def annotate(self,
                 insertions,
                 genes,
                 closest=False,
                 blacklist=None,
                 blacklist_col='gene_id'):
        """Annotates given insertions with predicted target genes."""

        matches = self.match(insertions, genes=genes)

        # Annotate and filter matches.
        matches = annotate_matches(matches, insertions=insertions, genes=genes)

        matches = filter_matches(
            matches,
            closest=closest,
            blacklist=blacklist,
            blacklist_col=blacklist_col)

        # Merge with insertion data.
        merged = pd.merge(insertions.values, matches, on='id', how='left')

        return InsertionSet(merged)

    def annotate_cis(self,
                     insertions,
                     genes,
                     cis_sites,
                     closest=False,
                     blacklist=None,
                     blacklist_col='gene_id'):
        """Annotates given insertions with predicted target genes,
           matched via CIS sites.
        """

        # TODO: Expand unstranded CIS sites?

        matches = self.match(cis_sites, genes=genes)

        # Annotate CIS matches with insertion ids.
        matches = pd.merge(
            matches.rename(columns={'id': 'cis_id'}),
            insertions.values['id', 'cis_id'],
            on='cis_id',
            how='left')

        # Annotate and filter matches.
        matches = annotate_matches(matches, insertions=insertions, genes=genes)

        matches = filter_matches(
            matches,
            closest=closest,
            blacklist=blacklist,
            blacklist_col=blacklist_col)

        # Merge with insertion data.
        merged = pd.merge(insertions.values, matches, on='id', how='left')

        return InsertionSet(merged)
