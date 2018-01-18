"""Module with functions for calling bowtie2."""

from pathlib import Path
import sys

import toolz

from . import util as shell


def bowtie2(read_paths,
            index_path,
            output_path,
            extra_options=None,
            verbose=False,
            logger=None):
    """
    Aligns reads to a reference genome using Bowtie2.

    Parameters
    ----------
    read_paths : List[Path]
        Path to input files containing reads.
    output_path : Path
        Output path for the aligned (and sorted) bam file.
    options : dict
        Dict of extra options to pass to Bowtie2. Should conform to the
        format expected by flatten_arguments.
    verbose : bool
        Whether to print output from bowtie2 to stderr.

    """

    assert isinstance(read_paths, tuple)
    assert 0 < len(read_paths) <= 2

    # Assemble options.
    input_options = {'-x': index_path}

    if len(read_paths) == 2:
        input_options['-1'] = read_paths[0]
        input_options['-2'] = read_paths[1]
    else:
        input_options['-U'] = read_paths[0]

    if any(ext in Path(read_paths[0]).suffixes for ext in {'.fa', '.fna'}):
        input_options['-f'] = True

    options = toolz.merge(extra_options or {}, input_options)

    # Build arguments.
    bowtie_args = ['bowtie2'] + shell.flatten_arguments(options)
    sort_args = ['samtools', 'sort', '-o', str(output_path), '-']

    # Run in piped fashion to avoid extra IO.
    processes = shell.run_piped([bowtie_args, sort_args])

    if verbose:
        # Print bowtie output to stderr for now.
        # TODO: Rewrite to use logging.
        print('', file=sys.stderr)
        stderr = processes[0].stderr.read().decode()
        print(stderr, file=sys.stderr)
