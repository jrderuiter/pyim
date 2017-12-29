"""Module containing the ShearSplink pipelines."""

import logging
from pathlib import Path

from cutadapt import seqio
import pandas as pd
import pysam

from pyim.external.cutadapt import cutadapt, cutadapt_summary
from pyim.external.bowtie2 import bowtie2
from pyim.external.util import flatten_arguments
from pyim.model import Insertion, InsertionSet
from pyim.util.path import shorten_path, extract_suffix

from .base import Pipeline, register_pipeline
from ..util import extract_insertions

DEFAULT_OVERLAP = 3
DEFAULT_ERROR_RATE = 0.1


class ShearSplinkPipeline(Pipeline):
    """ShearSplink pipeline.

    Analyzes (single-end) sequencing data that was prepared using the
    ShearSplink protocol. Sequence reads are expected to have the following
    structure::

        [Transposon][Genomic][Linker]

    Here, ``transposon`` refers to the flanking part of the transposon
    sequence, ``linker`` to the flanking linker sequence and ``genomic``
    to the genomic DNA located in between (which varies per insertion).
    The linker sequence is optional and may be omitted if the linker is not
    included in sequencing.

    The pipeline essentially performs the following steps:

        - If contaminants are provided, sequence reads are filtered
          (using Cutadapt) for the contaminant sequences.
        - The remaining reads are trimmed to remove the transposon and
          linker sequences, leaving only genomic sequences. Reads without
          the transposon/linker sequences are dropped, as we cannot be certain
          of their origin. (Note that the linker is optional and is only
          trimmed if a linker is given).
        - The genomic reads are aligned to the reference genome.
        - The resulting alignment is used to identify insertions.

    Note that this pipeline does **NOT** support multiplexed datasets (which is
    the default output of the ShearSplink protocol). For multiplexed datasets,
    use the ``MultiplexedShearSplinkPipeline``.

    Parameters
    ----------
    transposon_path : Path
        Path to the (flanking) transposon sequence (fasta).
    bowtie_index_path : Path
        Path to the bowtie index.
    linker_path : Path
        Path to the linker sequence (fasta).
    contaminant_path : Path
        Path to file containing contaminant sequences (fasta). If provided,
        sequences are filtered for these sequences before extracting genomic
        sequences for alignment.
    min_length : int
        Minimum length for genomic reads to be kept for alignment.
    min_support : int
        Minimum support for insertions to be kept in the final output.
    min_mapq : int
        Minimum mapping quality of alignments to be used for
        identifying insertions.
    merge_distance : int
        Maximum distance within which insertions are merged. Used to merge
        insertions that occur within close vicinity, which is typically due
        to slight variations in alignments.
    bowtie_options : Dict[str, Any]
        Dictionary of extra options for Bowtie.
    min_overlaps : Dict[str, int]
        Minimum overlap required to recognize the transposon, linker and
        contaminant sequences (see Cutadapts documentation for more
        information). Keys of the dictionary indicate to which sequence the
        overlap corresponds and should be one of the following: ``linker``,
        ``transposon`` or ``contaminant``.
    error_rates : Dict[str, float]
        Maximum error rate to use when recognizing transposon, linker and
        contaminant sequences (see Cutadapts documentation for more
        information). Keys should be the same as for ``min_overlaps``.

    """

    def __init__(self,
                 transposon_path,
                 bowtie_index_path,
                 linker_path=None,
                 contaminant_path=None,
                 min_length=15,
                 min_support=2,
                 min_mapq=23,
                 merge_distance=None,
                 bowtie_options=None,
                 min_overlaps=None,
                 error_rates=None):
        super().__init__()

        self._transposon_path = transposon_path
        self._linker_path = linker_path
        self._contaminant_path = contaminant_path

        self._index_path = bowtie_index_path

        self._min_length = min_length
        self._min_support = min_support
        self._min_mapq = min_mapq

        self._merge_distance = merge_distance
        self._bowtie_options = bowtie_options or {}

        self._min_overlaps = min_overlaps or {}
        self._error_rates = error_rates or {}

    @classmethod
    def configure_args(cls, parser):
        cls._setup_base_args(parser, paired=False)

        parser.description = 'ShearSplink pipeline'

        # Paths to various sequences.
        seq_options = parser.add_argument_group('Sequences')

        seq_options.add_argument(
            '--transposon',
            type=Path,
            required=True,
            help='Fasta file containing the transposon sequence.')

        seq_options.add_argument(
            '--contaminants',
            type=Path,
            default=None,
            help='Fasta file containing contaminant sequences.')

        seq_options.add_argument(
            '--linker',
            type=Path,
            default=None,
            help='Fasta file containing the linker sequence.')

        # Trimming options (used for cutadapt).
        trim_options = parser.add_argument_group('Trimming')

        trim_options.add_argument(
            '--min_length',
            type=int,
            default=15,
            help='Minimum length for (trimmed) genomic sequences.')

        trim_options.add_argument(
            '--contaminant_error',
            default=0.1,
            type=float,
            help='Maximum error rate for matching contaminants.')

        trim_options.add_argument(
            '--contaminant_overlap',
            default=3,
            type=int,
            help='Minimum overlap for matching contaminants.')

        trim_options.add_argument(
            '--transposon_error',
            default=0.1,
            type=float,
            help='Maximum error rate for matching the transposon.')

        trim_options.add_argument(
            '--transposon_overlap',
            default=3,
            type=int,
            help='Minimum overlap for matching the transposon.')

        trim_options.add_argument(
            '--linker_error',
            default=0.1,
            type=float,
            help='Maximum error rate for matching the linker.')

        trim_options.add_argument(
            '--linker_overlap',
            default=3,
            type=int,
            help='Minimum overlap for matching the linker.')

        align_options = parser.add_argument_group('Alignment')

        align_options.add_argument(
            '--bowtie_index',
            type=Path,
            required=True,
            help='Bowtie2 index to use for alignment.')

        align_options.add_argument(
            '--local',
            default=False,
            action='store_true',
            help='Use local alignment.')

        ins_options = parser.add_argument_group('Insertions')

        ins_options.add_argument(
            '--min_mapq',
            type=int,
            default=23,
            help=('Minimum mapping quality for reads '
                  'used to identify insertions.'))

        ins_options.add_argument(
            '--merge_distance',
            type=int,
            default=None,
            help=('Distance within which insertions (from same '
                  'sample) are merged.'))

        ins_options.add_argument(
            '--min_support',
            type=int,
            default=2,
            help='Minimum support for insertions.')

    @classmethod
    def _extract_args(cls, args):
        bowtie_options = {'--local': args.local}

        min_overlaps = {
            'contaminant': args.contaminant_overlap,
            'transposon': args.transposon_overlap,
            'linker': args.linker_overlap
        }

        error_rates = {
            'contaminant': args.contaminant_error,
            'transposon': args.transposon_error,
            'linker': args.linker_error
        }

        return dict(
            transposon_path=args.transposon,
            bowtie_index_path=args.bowtie_index,
            linker_path=args.linker,
            contaminant_path=args.contaminants,
            min_length=args.min_length,
            min_support=args.min_support,
            min_mapq=args.min_mapq,
            merge_distance=args.merge_distance,
            bowtie_options=bowtie_options,
            min_overlaps=min_overlaps,
            error_rates=error_rates)

    def run(self, read_path, output_dir, read2_path=None):
        if read2_path is not None:
            raise ValueError('Pipeline does not support paired-end data')

        logger = logging.getLogger()

        # Ensure output dir exists.
        output_dir.mkdir(exist_ok=True, parents=True)

        # Extract genomic sequences and align to reference.
        genomic_path = self._extract_genomic(read_path, output_dir, logger)
        alignment_path = self._align(genomic_path, output_dir, logger)

        # Extract insertions from bam file.
        bam_file = pysam.AlignmentFile(str(alignment_path))

        try:
            insertions = extract_insertions(
                iter(bam_file),
                func=_process_alignment,
                merge_dist=self._merge_distance,
                min_mapq=self._min_mapq,
                min_support=self._min_support,
                logger=logger)
        finally:
            bam_file.close()

        # Write insertions to output file.
        insertion_path = output_dir / 'insertions.txt'

        ins_frame = Insertion.to_frame(insertions)
        ins_frame.to_csv(str(insertion_path), sep='\t', index=False)

    def _extract_genomic(self, read_path, output_dir, logger):
        """Extracts the genomic part of sequence reads."""

        # Log parameters
        if logger is not None:
            logger.info('Extracting genomic sequences')
            logger.info('  %-18s: %s', 'Transposon',
                        shorten_path(self._transposon_path))
            logger.info('  %-18s: %s', 'Linker',
                        shorten_path(self._linker_path))
            logger.info('  %-18s: %s', 'Contaminants',
                        shorten_path(self._contaminant_path))
            logger.info('  %-18s: %s', 'Minimum length', self._min_length)

        # Get suffix to use for intermediate/genomic files.
        suffix = extract_suffix(read_path)

        # Track interim files for cleaning.
        interim_files = []

        if self._contaminant_path is not None:
            # Remove contaminants.
            contaminant_out_path = output_dir / (
                'trimmed_contaminant' + suffix)

            contaminant_opts = {
                '-g': 'file:' + str(self._contaminant_path),
                '--discard-trimmed': True,
                '-O': self._min_overlaps.get('contaminant', DEFAULT_OVERLAP),
                '-e': self._error_rates.get('contaminant', DEFAULT_ERROR_RATE)
            }

            process = cutadapt(read_path, contaminant_out_path,
                               contaminant_opts)

            if logger is not None:
                summary = cutadapt_summary(process.stdout, padding='   ')
                logger.info('Trimmed contaminant sequences' + summary)

            interim_files.append(contaminant_out_path)
        else:
            contaminant_out_path = read_path

        if self._linker_path is not None:
            # Remove linker.
            linker_out_path = output_dir / ('trimmed_linker' + suffix)
            linker_opts = {
                '-a': 'file:' + str(self._linker_path),
                '--discard-untrimmed': True,
                '-O': self._min_overlaps.get('linker', DEFAULT_OVERLAP),
                '-e': self._error_rates.get('linker', DEFAULT_ERROR_RATE)
            }

            process = cutadapt(contaminant_out_path, linker_out_path,
                               linker_opts)

            if logger is not None:
                summary = cutadapt_summary(process.stdout, padding='   ')
                logger.info('Trimmed linker sequence' + summary)

            interim_files.append(linker_out_path)
        else:
            linker_out_path = contaminant_out_path

        # Trim transposon and check minimum length.
        transposon_opts = {
            '-g': 'file:' + str(self._transposon_path),
            '--discard-untrimmed': True,
            '-O': self._min_overlaps.get('transposon', DEFAULT_OVERLAP),
            '-e': self._error_rates.get('transposon', DEFAULT_ERROR_RATE)
        }

        if self._min_length is not None:
            transposon_opts['--minimum-length'] = self._min_length

        genomic_path = output_dir / ('genomic' + suffix)
        process = cutadapt(linker_out_path, genomic_path, transposon_opts)

        if logger is not None:
            summary = cutadapt_summary(process.stdout, padding='   ')
            logger.info('Trimmed transposon sequence and filtered '
                        'for length' + summary)

        # Clean-up interim files.
        for file_path in interim_files:
            file_path.unlink()

        return genomic_path

    def _align(self, read_path, output_dir, logger):
        """Aligns genomic reads to the reference genome using Bowtie."""

        # Log parameters
        if logger is not None:
            logger.info('Aligning to reference')
            logger.info('  %-18s: %s', 'Reference',
                        shorten_path(self._index_path))
            logger.info('  %-18s: %s', 'Bowtie options',
                        flatten_arguments(self._bowtie_options))

        alignment_path = output_dir / 'alignment.bam'

        bowtie2(
            [read_path],
            index_path=self._index_path,
            output_path=alignment_path,
            options=self._bowtie_options,
            verbose=True)

        return alignment_path


register_pipeline(name='shearsplink', pipeline=ShearSplinkPipeline)


def _process_alignment(aln):
    """Analyzes an alignment to determine the tranposon/linker breakpoints."""
    ref = aln.reference_name

    if aln.is_reverse:
        transposon_pos = aln.reference_end
        linker_pos = aln.reference_start
        strand = -1
    else:
        transposon_pos = aln.reference_start
        linker_pos = aln.reference_end
        strand = 1

    return (ref, transposon_pos, strand), linker_pos


class MultiplexedShearSplinkPipeline(ShearSplinkPipeline):
    """ShearSplink pipeline supporting multiplexed reads.

    Analyzes multiplexed (single-end) sequencing data that was prepared using
    the ShearSplink protocol. Sequence reads are expected to have the following
    structure::

        [Barcode][Transposon][Genomic][Linker]

    Here, the ``transposon``, ``genomic`` and ``linker`` sequences are the
    same as for the ``ShearSplinkPipeline``. The ``barcode`` sequence is an
    index that indicates which sample the read originated for.

    Barcode sequences should be provided using the ``barcode_path`` argument.
    The optional ``barcode_mapping`` argument can be used to map barcodes to
    sample names.

    Parameters
    ----------
    transposon_path : Path
        Path to the (flanking) transposon sequence (fasta).
    bowtie_index_path : Path
        Path to the bowtie index.
    barcode_path :
        Path to barcode sequences (fasta).
    barcode_mapping : Path
        Path to a tsv file specifying a mapping from barcodes to sample names.
        Should contain ``sample`` and ``barcode`` columns.
    linker_path : Path
        Path to the linker sequence (fasta).
    contaminant_path : Path
        Path to file containing contamintant sequences (fasta). If provided,
        sequences are filtered for these sequences before extracting genomic
        sequences for alignment.
    min_length : int
        Minimum length for genomic reads to be kept for alignment.
    min_support : int
        Minimum support for insertions to be kept in the final output.
    min_mapq : int
        Minimum mapping quality of alignments to be used for
        identifying insertions.
    merge_distance : int
        Maximum distance within which insertions are merged. Used to merge
        insertions that occur within close vicinity, which is typically due
        to slight variations in alignments.
    bowtie_options : Dict[str, Any]
        Dictionary of extra options for Bowtie.
    min_overlaps : Dict[str, int]
        Minimum overlap required to recognize the transposon, linker and
        contamintant sequences (see Cutadapts documentation for more
        information). Keys of the dictionary indicate to which sequence the
        overlap corresponds and should be one of the following: ``linker``,
        ``transposon`` or ``contaminant``.
    error_rates : Dict[str, float]
        Maximum error rate to use when recognizing transposon, linker and
        contamintant sequences (see Cutadapts documentation for more
        information). Keys should be the same as for ``min_overlaps``.

    """

    def __init__(self,
                 transposon_path,
                 bowtie_index_path,
                 barcode_path,
                 barcode_mapping=None,
                 linker_path=None,
                 contaminant_path=None,
                 min_length=15,
                 min_support=2,
                 min_mapq=23,
                 merge_distance=0,
                 bowtie_options=None,
                 min_overlaps=None,
                 error_rates=None):
        super().__init__(
            transposon_path=transposon_path,
            bowtie_index_path=bowtie_index_path,
            linker_path=linker_path,
            contaminant_path=contaminant_path,
            min_length=min_length,
            min_support=min_support,
            min_mapq=min_mapq,
            merge_distance=merge_distance,
            bowtie_options=bowtie_options,
            min_overlaps=min_overlaps,
            error_rates=error_rates)

        self._barcode_path = barcode_path
        self._barcode_mapping = barcode_mapping

    @classmethod
    def configure_args(cls, parser):
        super().configure_args(parser)

        parser.add_argument('--barcodes', required=True, type=Path)
        parser.add_argument(
            '--barcode_mapping', required=False, type=Path, default=None)

    @classmethod
    def _extract_args(cls, args):
        arg_dict = super()._extract_args(args)

        if args.barcode_mapping is not None:
            map_df = pd.read_csv(args.barcode_mapping, sep='\t')
            arg_dict['barcode_mapping'] = dict(
                zip(map_df['barcode'], map_df['sample']))
        else:
            arg_dict['barcode_mapping'] = None

        arg_dict['barcode_path'] = args.barcodes

        return arg_dict

    def run(self, read_path, output_dir, read2_path=None):
        if read2_path is not None:
            raise ValueError('Pipeline does not support paired-end data')

        logger = logging.getLogger()

        # Ensure output dir exists.
        output_dir.mkdir(exist_ok=True, parents=True)

        # Extract genomic sequences and align to reference.
        genomic_path = self._extract_genomic(read_path, output_dir, logger)
        alignment_path = self._align(genomic_path, output_dir, logger)

        # Map reads to specific barcodes/samples.
        logger.info('Extracting barcode/sample mapping')
        logger.info('  %-18s: %s', 'Barcodes',
                    shorten_path(self._barcode_path))
        read_map = self._get_barcode_mapping(read_path)

        # Extract insertions from bam file.
        bam_file = pysam.AlignmentFile(str(alignment_path))

        try:
            insertions = extract_insertions(
                iter(bam_file),
                func=_process_alignment,
                group_func=lambda aln: read_map.get(aln.query_name, None),
                merge_dist=self._merge_distance,
                min_mapq=self._min_mapq,
                min_support=self._min_support,
                logger=logger)
        finally:
            bam_file.close()

        # Write insertions to output file.
        insertion_path = output_dir / 'insertions.txt'

        insertion_set = InsertionSet.from_tuples(insertions)
        insertion_set.to_csv(str(insertion_path), sep='\t', index=False)

    def _get_barcode_mapping(self, read_path):
        # Read barcode sequences.
        with seqio.open(str(self._barcode_path)) as barcode_file:
            barcodes = list(barcode_file)

        # Extract read --> barcode mapping.
        with seqio.open(str(read_path)) as reads:
            return _extract_barcode_mapping(reads, barcodes,
                                            self._barcode_mapping)


register_pipeline(
    name='shearsplink-multiplexed', pipeline=MultiplexedShearSplinkPipeline)


def _extract_barcode_mapping(reads, barcodes, barcode_mapping=None):

    # Create barcode/sample dict.
    barcode_dict = {bc.name: bc.sequence for bc in barcodes}

    if barcode_mapping is not None:
        barcode_dict = {
            sample: barcode_dict[barcode]
            for barcode, sample in barcode_mapping.items()
        }

    # Build mapping.
    mapping = {}

    for read in reads:
        # Check each barcode for match in read.
        matched = [k for k, v in barcode_dict.items() if v in read.sequence]

        if len(matched) == 1:
            # Record single matches.
            name = read.name.split()[0]
            mapping[name] = matched[0]
        elif len(matched) > 1:
            logging.warning('Skipping %s due to multiple matching barcodes',
                            read.name.split()[0])

    return mapping
