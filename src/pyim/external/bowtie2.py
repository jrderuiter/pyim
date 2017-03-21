"""Module with functions for calling bowtie2."""

import sys

from . import util as shell


def bowtie2(in1_paths,
            index_path,
            output_path,
            options=None,
            in2_paths=None,
            verbose=False):
    """
    Aligns reads to a reference genome using Bowtie2.

    Parameters
    ----------
    in1_paths : List[Path]
        Path to input files containings reads. For single read data,
        a list of Paths is expected. For paired-end sequencing data,
        Paths should be passed as a tuple of lists, in which the first
        element is taken as #1 mates and the second as #2 mates.
    output_path : Path
        Output path for the aligned (and sorted) bam file.
    options : dict
        Dict of extra options to pass to Bowtie2.
    """

    # Ensure we have a copy of options to work on.
    options = dict(options) if options is not None else {}

    # Inject inputs + index into options.
    if in2_paths is not None:
        options['-1'] = ','.join(str(fp) for fp in in1_paths)
        options['-2'] = ','.join(str(fp) for fp in in2_paths)
    else:
        options['-U'] = ','.join(str(fp) for fp in in1_paths)

    if any(ext in in1_paths[0].suffixes for ext in {'.fa', '.fna'}):
        options['-f'] = True

    options['-x'] = str(index_path)

    # Build bowtie2 arguments.
    bowtie_args = ['bowtie2'] + shell.flatten_arguments(options)

    # Sort arguments for samtools.
    sort_args = ['samtools', 'sort', '-o', str(output_path), '-']

    # Run in piped fashion to avoid extra IO.
    processes = shell.run_piped([bowtie_args, sort_args])

    if verbose:
        # Print bowtie output to stderr for now.
        # TODO: Rewrite to use logging.
        print('', file=sys.stderr)
        stderr = processes[0].stderr.read().decode()
        print(stderr, file=sys.stderr)
