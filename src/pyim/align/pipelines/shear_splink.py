"""Module containing the ShearSplink pipelines."""

import logging
from pathlib import Path

from cutadapt import seqio
import pandas as pd
import pysam

from pyim.external.cutadapt import cutadapt, cutadapt_summary
from pyim.external.bowtie2 import bowtie2
from pyim.external.util import flatten_arguments
from pyim.model import Insertion
from pyim.util.path import shorten_path, extract_suffix

from .base import Pipeline, register_pipeline
from ..util import extract_insertions

DEFAULT_OVERLAP = 3
DEFAULT_ERROR_RATE = 0.1


class ShearSplinkPipeline(Pipeline):
    """ShearSplink pipeline."""

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

        parser.add_argument('--transposon', type=Path, required=True)
        parser.add_argument('--bowtie_index', type=Path, required=True)

        parser.add_argument('--contaminants', type=Path, default=None)
        parser.add_argument('--linker', type=Path, default=None)

        parser.add_argument('--min_length', type=int, default=15)
        parser.add_argument('--min_support', type=int, default=2)
        parser.add_argument('--min_mapq', type=int, default=23)
        parser.add_argument('--merge_distance', type=int, default=None)

        parser.add_argument('--local', default=False, action='store_true')

        parser.add_argument('--contaminant_error', default=0.1, type=float)
        parser.add_argument('--transposon_error', default=0.1, type=float)
        parser.add_argument('--linker_error', default=0.1, type=float)

        parser.add_argument('--contaminant_overlap', default=3, type=int)
        parser.add_argument('--transposon_overlap', default=3, type=int)
        parser.add_argument('--linker_overlap', default=3, type=int)

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
    """ShearSplink pipeline with multiplexed reads."""

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

        ins_frame = Insertion.to_frame(insertions)
        ins_frame.to_csv(str(insertion_path), sep='\t', index=False)

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
        barcode_dict = {sample: barcode_dict[barcode]
                        for barcode, sample in barcode_mapping.items()}

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
