"""Module providing base functionality for insertion identification
pipelines."""

import abc
from pathlib import Path

from pyim.main import Command


class Aligner(abc.ABC):
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

    @abc.abstractmethod
    def run(self, read_paths, work_dir=None):
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


class AlignerCommand(Command):
    """Base aligner command."""
    pass


class SingleEndCommand(AlignerCommand):
    """Base command class for single-end aligners."""

    def configure(self, parser):
        parser.add_argument('--reads', required=True, type=Path)
        parser.add_argument('--output', required=True, type=Path)


class PairedEndCommand(AlignerCommand):
    """Base command class for paired-end aligners."""

    def configure(self, parser):
        parser.add_argument('--reads', nargs=2, required=True)
        parser.add_argument('--output', required=True, type=Path)
