from abc import ABC, abstractclassmethod, abstractmethod
from pathlib import Path

_registry = {}


def register_caller(name, cis_caller):
    _registry[name] = cis_caller


def get_callers():
    return dict(_registry)


class CisCaller(ABC):
    def __init__(self):
        pass

    @abstractclassmethod
    def configure_args(cls, parser):
        parser.add_argument('--insertions', type=Path, required=True)
        parser.add_argument('--output', type=Path, required=True)
        parser.add_argument(
            '--output_sites', type=Path, required=False, default=None)

    @abstractclassmethod
    def from_args(cls, args):
        raise NotImplementedError()

    @abstractmethod
    def call(self, insertions):
        """Calls CIS sites for insertions.

        Parameters:
            insertions (iterable[Insertion])

        Returns:
            iterable[Insertions], iterable[CisSites]

        """

        raise NotImplementedError()