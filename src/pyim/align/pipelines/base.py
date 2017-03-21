import abc
from pathlib import Path

_registry = {}


def register_pipeline(name, pipeline):
    _registry[name] = pipeline


def get_pipelines():
    return dict(_registry)


class Pipeline(abc.ABC):
    """Base pipeline class."""

    def __init__(self):
        pass

    @abc.abstractclassmethod
    def configure_args(cls, parser):
        """Configures argument parser for the pipeline."""
        parser.add_argument('--reads', type=Path, required=True)
        parser.add_argument('--reads2', type=Path, required=False)
        parser.add_argument('--output_dir', type=Path, required=True)

    @classmethod
    def from_args(cls, args):
        """Builds a pipeline instance from the given arguments."""
        return cls(**cls._extract_args(args))

    @abc.abstractclassmethod
    def _extract_args(cls, args):
        """Extract arguments from args for from_args."""

    @abc.abstractmethod
    def run(self, reads_path, output_dir, reads2_path=None):
        """Runs the pipeline with the given input."""
