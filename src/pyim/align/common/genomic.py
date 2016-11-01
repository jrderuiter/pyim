import logging

from pyim.util.path import build_path

from pyim.external.cutadapt import cutadapt, cutadapt_summary

DEFAULT_OVERLAP = 3
DEFAULT_ERROR_RATE = 0.1


def extract_genomic(reads_path,
                    output_path,
                    transposon_path,
                    linker_path=None,
                    contaminant_path=None,
                    min_length=None,
                    min_overlaps=None,
                    error_rates=None):
    """Extracts genomic sequences from single-read data."""

    logger = logging.getLogger()

    min_overlaps = min_overlaps or {}
    error_rates = error_rates or {}

    # Ensure output dir exists.
    output_path.parent.mkdir(exist_ok=True)

    # Track interim files for cleaning.
    interim_files = []

    if contaminant_path is not None:
        # Remove contaminants.
        contaminant_out_path = build_path(output_path, suffix='.contaminant')
        contaminant_opts = {
            '-g': 'file:' + str(contaminant_path),
            '--discard-trimmed': True,
            '-O': min_overlaps.get('contaminant', DEFAULT_OVERLAP),
            '-e': error_rates.get('contaminant', DEFAULT_ERROR_RATE)
        }

        p = cutadapt(reads_path, contaminant_out_path, contaminant_opts)
        logger.info('Trimmed contaminant sequences' +
                    cutadapt_summary(p.stdout)) # yapf: disable

        interim_files.append(contaminant_out_path)
    else:
        contaminant_out_path = reads_path

    if linker_path is not None:
        # Remove linker.
        linker_out_path = build_path(output_path, suffix='.linker')
        linker_opts = {
            '-a': 'file:' + str(linker_path),
            '--discard-untrimmed': True,
            '-O': min_overlaps.get('linker', DEFAULT_OVERLAP),
            '-e': error_rates.get('linker', DEFAULT_ERROR_RATE)
        }

        p = cutadapt(contaminant_out_path, linker_out_path, linker_opts)
        logger.info('Trimmed linker sequence' +
                    cutadapt_summary(p.stdout)) # yapf: disable

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

    p = cutadapt(linker_out_path, output_path, transposon_opts)
    logger.info('Trimmed transposon sequence and filtered for length' +
                cutadapt_summary(p.stdout)) # yapf: disable

    # Clean-up interim files.
    for fp in interim_files:
        fp.unlink()


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
                 in2_path=in2_path, out2_path=out2_path) # yapf: disable

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
                 in2_path=cont2_path, out2_path=link2_path) # yapf: disable

        interim_files += [link1_path, link2_path]
    else:
        link1_path, link2_path = cont1_path, cont2_path

    # Trim transposon and check minimum length.
    transposon_opts = {'-g': 'file:' + str(transposon_path),
                       '--discard-untrimmed': True}

    if min_length is not None:
        transposon_opts['--minimum-length'] = min_length

    cutadapt(link1_path, out1_path, transposon_opts,
             in2_path=link2_path, out2_path=out2_path) # yapf: disable

    # Clean-up intermediary files.
    for fp in interim_files:
        fp.unlink()
