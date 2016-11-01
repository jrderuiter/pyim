from abc import ABC, abstractclassmethod, abstractmethod, abstractproperty
import functools
import itertools
import operator
from pathlib import Path

from frozendict import frozendict
import numpy as np
import toolz

from pyim.model import Insertion, CisSite
from ..metadata import add_metadata

_registry = {}


def register_annotator(name, annotator):
    _registry[name] = annotator


def get_annotators():
    return dict(_registry)


class Annotator(ABC):
    def __init__(self):
        pass

    @classmethod
    def configure_args(cls, parser):
        parser.add_argument('--insertions', type=Path, required=True)
        parser.add_argument('--output', type=Path, required=True)

    @classmethod
    def from_args(cls, args):
        return cls(**cls.parse_args(args))

    @abstractclassmethod
    def parse_args(cls, args):
        raise NotImplementedError()

    @abstractmethod
    def annotate(self, insertions):
        raise NotImplementedError()

    @abstractproperty
    def gtf(self):
        raise NotImplementedError()


class CisAnnotator(Annotator):
    def __init__(self, *args, cis_sites=None, drop_cis_id=False, **kwargs):
        super().__init__(*args, **kwargs)

        self._cis_sites = self._preprocess_sites(cis_sites)
        self._drop_cis_id = drop_cis_id

    def _preprocess_sites(self, cis_sites):
        """Pre-process cis sites, fixing unstrandedness etc."""

        # Copy CISs that are unstranded to both strands.
        return list(self._expand_unstranded_sites(cis_sites))

    @staticmethod
    def _expand_unstranded_sites(cis_sites):
        for cis in cis_sites:
            if np.isnan(cis.strand) or cis.strand is None:
                yield cis._replace(strand=1)
                yield cis._replace(strand=-1)
            else:
                yield cis

    @classmethod
    def configure_args(cls, parser):
        super().configure_args(parser)
        parser.add_argument('--cis_sites', default=None, type=Path)
        parser.add_argument(
            '--drop_cis_id', default=False, action='store_true')

    @classmethod
    def parse_args(cls, args):
        parsed = super().parse_args(args)

        if args.cis_sites is not None:
            cis_cols = ['id', 'chromosome', 'position', 'strand']
            cis_sites = list(
                CisSite.from_csv(
                    args.cis_sites, usecols=cis_cols, sep='\t'))
        else:
            cis_sites = None

        cis_args = {'cis_sites': cis_sites, 'drop_cis_id': args.drop_cis_id}

        return toolz.merge(parsed, cis_args)

    def annotate(self, insertions):
        if self._cis_sites is None:
            yield from super().annotate(insertions)
        else:
            yield from self._annotate_cis(insertions)

    def _annotate_cis(self, insertions):
        # Annotate CIS sites.
        annotated_sites = super().annotate(self._cis_sites)

        # Create CIS --> gene map using annotations.
        id_getter = operator.attrgetter('id')
        annotated_sites = sorted(annotated_sites, key=id_getter)
        grouped_sites = itertools.groupby(annotated_sites, key=id_getter)

        cis_map = {
            cis_id: {(item.metadata['gene_name'], item.metadata['gene_id'])
                     for item in group if 'gene_name' in item.metadata}
            for cis_id, group in grouped_sites
        }

        # Annotate insertions, drop any duplicates and add metadata.
        annotated_ins = set(self._annotate_insertions(insertions, cis_map))
        annotated_ins = add_metadata(annotated_ins, reference_gtf=self.gtf)

        yield from annotated_ins

    def _annotate_insertions(self, insertions, cis_map):
        for insertion in insertions:
            genes = cis_map.get(insertion.metadata['cis_id'], set())

            if len(genes) > 0:
                for gene_name, gene_id in genes:
                    metadata = {'gene_id': gene_id, 'gene_name': gene_name}
                    metadata = toolz.merge(insertion.metadata, metadata)

                    if self._drop_cis_id:
                        metadata.pop('cis_id')

                    yield insertion._replace(metadata=frozendict(metadata))
            else:
                if self._drop_cis_id:
                    metadata = dict(insertion.metadata)
                    metadata.pop('cis_id')
                    yield insertion._replace(metadata=frozendict(metadata))
                else:
                    yield insertion

    @property
    def gtf(self):
        return super().gtf
