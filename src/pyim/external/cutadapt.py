"""Module with functions for calling cutadapt."""

import itertools
from pathlib import Path
import shutil

import pyfaidx

from . import util as shell


def cutadapt(in1_path, out1_path, options, in2_path=None, out2_path=None):
    """Runs cutadapt using the given options."""

    in1_path = in1_path or '-'
    options = dict(options) if options is not None else {}

    if out1_path is not None:
        options['-o'] = str(out1_path)

    if out2_path is not None:
        options['-p'] = str(out2_path)

    cmdline_args = shell.flatten_arguments(options)
    cmdline_args = ['cutadapt'] + cmdline_args + [str(in1_path)]

    if in2_path is not None:
        cmdline_args += [str(in2_path)]

    return shell.run(cmdline_args)


def demultiplex_samples(reads_path,
                        output_dir,
                        barcode_path,
                        error_rate=0.0,
                        sample_mapping=None):
    """
    De-multiplexes reads into separate sample/barcode files.

    Parameters
    ----------
    reads_path : Path
        Path to the input reads file (in fasta/fastq format).
    output_dir : Path
        Output directory to which the de-multiplexed files will be written.
    barcode_path : Path
        Path to fasta file containing barcode sequences.
    sample_mapping : dict
        Dict mapping barcodes to samples.

    Returns
    -------
    dict[Path]
        Returns dict mapping samples to the respective demultiplexed file.
    """

    if sample_mapping is None:
        # Directly de-multiplex using barcodes.
        sample_paths = _demultiplex(
            reads_path, output_dir, barcode_path, error_rate=error_rate)
    else:
        # First demultiplex to barcodes in temp dir.
        tmp_dir = output_dir / '_barcodes'
        barcode_paths = _demultiplex(
            reads_path, tmp_dir, barcode_path, error_rate=error_rate)

        # Then rename files using mapping and delete files for unused barcodes.
        sample_paths = {}

        for barcode, sample in sample_mapping.items():
            barcode_path = barcode_paths[barcode]
            sample_path = output_dir / (sample + barcode_path.suffixes[-1])

            if barcode_path.exists():
                shutil.move(str(barcode_path), str(sample_path))
            else:
                # Create empty output if nothing was extracted for barcode.
                sample_path.touch()

            sample_paths[sample] = sample_path

        shutil.rmtree(str(tmp_dir))

    return sample_paths


def _demultiplex(reads_path, output_dir, barcode_path, error_rate):
    """Runs cutadapt to de-multiplex reads into seperate files per barcode."""

    output_dir.mkdir(parents=True)

    # De-multiplex using cutadapt.
    options = {'-g': 'file:' + str(barcode_path),
               '--discard-untrimmed': True,
               '-e': error_rate}
    output_base = output_dir / ('{name}' + reads_path.suffixes[-1])
    cutadapt(reads_path, output_base, options=options)

    # Identify output files.
    barcode_keys = pyfaidx.Fasta(str(barcode_path)).keys()
    output_paths = {bc: Path(str(output_base).format(name=bc))
                    for bc in barcode_keys}

    return output_paths


def cutadapt_summary(stdstream):
    sections = _split_log_sections(stdstream.read().decode())
    return '\n'.join([' '] + sections['=== Summary ==='])


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
