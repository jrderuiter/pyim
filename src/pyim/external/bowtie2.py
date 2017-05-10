"""Module with functions for calling bowtie2."""

import sys

from . import util as shell


def bowtie2(read_paths,
            index_path,
            output_path,
            options=None,
            read2_paths=None,
            verbose=False):
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
    read2_paths : List[Path]
        Path to input files containing the second end (for paired-end data).
    verbose : bool
        Whether to print output from bowtie2 to stderr.

    """

    # Ensure we have a copy of options to work on.
    options = dict(options) if options is not None else {}

    # Inject inputs + index into options.
    if read2_paths is not None:
        options['-1'] = ','.join(str(fp) for fp in read_paths)
        options['-2'] = ','.join(str(fp) for fp in read2_paths)
    else:
        options['-U'] = ','.join(str(fp) for fp in read_paths)

    if any(ext in read_paths[0].suffixes for ext in {'.fa', '.fna'}):
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
