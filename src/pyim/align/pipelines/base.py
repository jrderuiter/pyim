"""Module providing base functionality for insertion identification
pipelines."""

import abc
from pathlib import Path

_registry = {}


def register_pipeline(name, pipeline):
    """Registers a pipeline class under the given name.

    Parameters
    ----------
    name : str
        Name to use for the pipeline.
    pipeline : Pipeline
        The pipeline class.

    """
    _registry[name] = pipeline


def get_pipelines():
    """Returns a dict of the available pipelines, indexed by pipeline name.

    Returns
    -------
    Dict[str, Pipeline]
        Available pipelines.

    """
    return dict(_registry)


class Pipeline(abc.ABC):
    """Base pipeline class.

    Pipeline classes implement analyses that derive transposon insertion sites
    from sequencing data obtained by targeted (DNA) sequencing of the insertion
    sites.

    The main interface of the class is the ``run`` method, whose main
    arguments are the paths to the sequence read files and the output
    directory. After completion, the output directory contains an
    ``insertions.txt`` output file, describing the location of the identified
    insertion sites, and any optional extra intermediate/output files.

    Each pipeline should also provide implementations for the
    ``configure_args`` and ``from_args`` methods, which are used to instantiate
    pipelines from command line arguments as part of the ``pyim-align`` command.

    """

    def __init__(self):
        pass

    @abc.abstractclassmethod
    def configure_args(cls, parser):
        """Configures argument parser for the pipeline.

        Parameters
        ----------
        parser : ArgumentParser
            ArgumentParser to configure.

        """

    @classmethod
    def _setup_base_args(cls, parser, paired=False):
        parser.add_argument('--reads', type=Path, required=True)

        if paired:
            parser.add_argument('--reads2', type=Path, required=False)

        parser.add_argument('--output_dir', type=Path, required=True)

    @classmethod
    def from_args(cls, args):
        """Builds a pipeline instance from the given arguments.

        Parameters
        ----------
        args : Namespace
            Parsed arguments from argparser.

        Returns
        -------
        Pipeline
            Instantiated pipeline instance.

        """
        return cls(**cls._extract_args(args))

    @abc.abstractclassmethod
    def _extract_args(cls, args):
        """Extracts arguments from args for from_args.

        Returns arguments as a dict of Dict[str, Any], specifying values for
        the various parameters of the corresponding pipeline class.

        Parameters
        ----------
        args : Namespace
            Parsed arguments from argparser.

        Returns
        -------
        Dict[str, Any]
            Dictionary of pipeline parameters.

        """

    @abc.abstractmethod
    def run(self, read_path, output_dir, read2_path=None):
        """Runs the pipeline, producing a table of identified insertions.

        Parameters
        ----------
        read_path : Path
            Path to sequence reads. For paired-end data, this should refer
            to the first read of the mate pair.
        output_dir : Path
            Path to the output directory.
        read2_path : Path
            Optional path to the second read of the mate pair (for paired-end)
            sequencing data. Only used in pipelines that support paired-end
            sequencing data.

        """
