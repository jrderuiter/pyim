"""Module providing base functionality for insertion identification
pipelines."""

import abc
from pathlib import Path


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
    def extract_insertions(self, read_paths, output_prefix=None):
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
