"""Module with functions for calling cutadapt."""

import itertools
from pathlib import Path

from . import util as shell


def cutadapt(input_paths,
             output_paths,
             options=None,
             threads=1,
             verbose=False,
             logger=None):
    """Runs cutadapt using the given options."""

    assert isinstance(input_paths, tuple)
    assert isinstance(output_paths, tuple)
    assert len(input_paths) == len(output_paths)

    # if (read2_path is None) != (out2_path is None):
    #     raise ValueError('Both read2_path and out2_path must be specified '
    #                      'for paired-end sequencing data')

    # Ensure output path exists.
    output_dir = Path(output_paths[0]).parent
    output_dir.mkdir(exist_ok=True, parents=True)

    # Assemble options.
    options = dict(options)

    options['-o'] = str(output_paths[0])
    options['-j'] = threads

    if len(output_paths) == 2:
        options['-p'] = str(output_paths[1])

    # Assemble arg list.
    inputs = [str(fp) for fp in input_paths]
    options_flat = shell.flatten_arguments(options)
    cmdline_args = ['cutadapt'] + options_flat + inputs

    process = shell.run(cmdline_args)

    # Log result if logger given.
    if verbose:
        summary = extract_summary_from_log(process.stdout)

        if logger is not None:
            logger.info(summary)
        else:
            print(summary)


def extract_summary_from_log(stdstream, padding=''):
    """Extract trim summary from log."""

    sections = _split_log_sections(stdstream.read().decode())
    delim = '\n' + padding

    return padding + delim.join([''] + sections['=== Summary ==='])


def _split_log_sections(log_str):
    return dict(_iter_log_sections(log_str.split('\n')))


def _iter_log_sections(lines):
    grouped = itertools.groupby(lines, lambda line: line.startswith('==='))
    group_iter = (x[1] for x in grouped)

    yield 'Header', list(next(group_iter))

    for name in group_iter:
        header = next(name).strip()
        lines = list(next(group_iter))
        yield header, lines


def _parse_summary_section(lines):
    lines = (line.strip() for line in lines)

    stats = {}
    for line in lines:
        if line:
            key, value = line.split(':')

            value = value.strip().split()[0]
            value = value.replace(',', '')

            stats[key] = int(value)

    return stats
