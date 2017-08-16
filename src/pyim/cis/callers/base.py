 import abc
from pathlib import Path

import toolz

from pyim.main import Command
from pyim.model import Insertion, CisSite
from pyim.vendor.frozendict import frozendict


class CisCaller(abc.ABC):
    """Base CisCaller class."""

    def __init__(self):
        pass

    @abc.abstractmethod
    def run(self, insertions):
        """Calls CIS sites for insertions.

        Parameters:
            insertions (iterable[Insertion])

        Returns:
            iterable[Insertions], iterable[CisSites]

        """

        raise NotImplementedError()


class CisCallerCommand(Command):
    """Base CisCaller command."""

    def configure(self, parser):
        parser.add_argument('--insertions', type=Path, required=True)
        parser.add_argument('--output', type=Path, required=True)
        parser.add_argument(
            '--output_sites', type=Path, required=False, default=None)

    @staticmethod
    def _annotate_insertions(insertions, cis_mapping):
        """Annotates insertions with CIS sites using given mapping."""

        for insertion in insertions:
            cis_ids = cis_mapping.get(insertion.id, set())

            for cis_id in cis_ids:
                cis_metadata = {'cis_id': cis_id}
                new_metadata = toolz.merge(insertion.metadata, cis_metadata)
                yield insertion._replace(metadata=frozendict(new_metadata))

    @staticmethod
    def _write_outputs(insertions, cis_sites, args):
        Insertion.to_csv(args.output, insertions, sep='\t', index=False)

        if args.output_sites is None:
            cis_path = args.output.with_suffix('.sites.txt')
        else:
            cis_path = args.output_sites

        CisSite.to_csv(cis_path, cis_sites, sep='\t', index=False)
