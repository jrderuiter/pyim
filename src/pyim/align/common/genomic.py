import logging

from pyim.util.path import build_path

from pyim.external.cutadapt import cutadapt, cutadapt_summary

DEFAULT_OVERLAP = 3
DEFAULT_ERROR_RATE = 0.1


def extract_genomic(reads_path,
                    output_dir,
                    transposon_path,
                    linker_path=None,
                    contaminant_path=None,
                    min_length=None,
                    min_overlaps=None,
                    error_rates=None):
    """Extracts genomic sequences from single-read data.

    Process reads of the following structure:

        [Transposon-Genomic-Linker]

    """

    logger = logging.getLogger()

    min_overlaps = min_overlaps or {}
    error_rates = error_rates or {}

    # Ensure output dir exists.
    output_dir.mkdir(exist_ok=True)

    suffix = _extract_suffix(reads_path)

    # Track interim files for cleaning.
    interim_files = []

    if contaminant_path is not None:
        # Remove contaminants.
        contaminant_out_path = output_dir / ('filt_contaminant' + suffix)
        contaminant_opts = {
            '-g': 'file:' + str(contaminant_path),
            '--discard-trimmed': True,
            '-O': min_overlaps.get('contaminant', DEFAULT_OVERLAP),
            '-e': error_rates.get('contaminant', DEFAULT_ERROR_RATE)
        }

        process = cutadapt(reads_path, contaminant_out_path, contaminant_opts)
        logger.info('Trimmed contaminant sequences' +
                    cutadapt_summary(process.stdout)) # yapf: disable

        interim_files.append(contaminant_out_path)
    else:
        contaminant_out_path = reads_path

    if linker_path is not None:
        # Remove linker.
        linker_out_path = output_dir / ('filt_linker' + suffix)
        linker_opts = {
            '-a': 'file:' + str(linker_path),
            '--discard-untrimmed': True,
            '-O': min_overlaps.get('linker', DEFAULT_OVERLAP),
            '-e': error_rates.get('linker', DEFAULT_ERROR_RATE)
        }

        process = cutadapt(contaminant_out_path, linker_out_path, linker_opts)
        logger.info('Trimmed linker sequence' +
                    cutadapt_summary(process.stdout)) # yapf: disable

        interim_files.append(linker_out_path)
    else:
        linker_out_path = contaminant_out_path

    # Trim transposon and check minimum length.
    transposon_opts = {
        '-g': 'file:' + str(transposon_path),
        '--discard-untrimmed': True,
        '-O': min_overlaps.get('transposon', DEFAULT_OVERLAP),
        '-e': error_rates.get('transposon', DEFAULT_ERROR_RATE)
    }

    if min_length is not None:
        transposon_opts['--minimum-length'] = min_length

    genomic_path = output_dir / ('genomic' + suffix)
    process = cutadapt(linker_out_path, genomic_path, transposon_opts)
    logger.info('Trimmed transposon sequence and filtered for length' +
                cutadapt_summary(process.stdout)) # yapf: disable

    # Clean-up interim files.
    for file_path in interim_files:
        file_path.unlink()

    return genomic_path


def _extract_suffix(file_path):
    if file_path.suffixes[-1] == '.gz':
        suffix = ''.join(file_path.suffixes[-2:])
    else:
        suffix = file_path.suffixes[-1]
    return suffix


def extract_genomic_paired(reads_paths,
                           output_paths,
                           transposon_path,
                           linker_path=None,
                           contaminant_path=None,
                           min_length=None):
    """Extracts genomic sequences from paired-end data."""

    # Extract file paths.
    in1_path, in2_path = reads_paths
    out1_path, out2_path = output_paths

    # Ensure output dirs exists.
    out1_path.parent.mkdir(exist_ok=True)
    out2_path.parent.mkdir(exist_ok=True)

    # Track interim files.
    interim_files = []

    if contaminant_path is not None:
        # Remove contaminants.
        cont1_path = build_path(out1_path, suffix='.contaminant')
        cont2_path = build_path(out2_path, suffix='.contaminant')

        contaminant_opts = {'-g': 'file:' + str(contaminant_path),
                            '--discard-trimmed': True}
        cutadapt(in1_path, cont1_path, contaminant_opts,
                 reads2_path=in2_path, out2_path=out2_path) # yapf: disable

        interim_files += [cont1_path, cont2_path]
    else:
        cont1_path, cont2_path = in1_path, in2_path

    if linker_path is not None:
        # Remove linker.
        link1_path = build_path(out1_path, suffix='.linker')
        link2_path = build_path(out2_path, suffix='.linker')

        linker_opts = {'-A': 'file:' + str(linker_path),
                       '--discard-untrimmed': True}
        cutadapt(cont1_path, link1_path, linker_opts,
                 reads2_path=cont2_path, out2_path=link2_path) # yapf: disable

        interim_files += [link1_path, link2_path]
    else:
        link1_path, link2_path = cont1_path, cont2_path

    # Trim transposon and check minimum length.
    transposon_opts = {'-g': 'file:' + str(transposon_path),
                       '--discard-untrimmed': True}

    if min_length is not None:
        transposon_opts['--minimum-length'] = min_length

    cutadapt(link1_path, out1_path, transposon_opts,
             reads2_path=link2_path, out2_path=out2_path) # yapf: disable

    # Clean-up intermediary files.
    for fp in interim_files:
        fp.unlink()
