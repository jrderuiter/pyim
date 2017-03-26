"""Module containing the nextera pipeline."""

import logging
from pathlib import Path

import pysam
import toolz

from pyim.external.bowtie2 import bowtie2
from pyim.external.cutadapt import cutadapt, cutadapt_summary
from pyim.external.util import flatten_arguments
from pyim.model import Insertion
from pyim.util.path import shorten_path, extract_suffix

from .base import Pipeline, register_pipeline
from ..util import extract_insertions


class NexteraPipeline(Pipeline):
    """Nextera-based transposon pipeline."""

    def __init__(self,
                 transposon_path,
                 bowtie_index_path,
                 bowtie_options=None,
                 min_length=15,
                 min_support=2,
                 min_mapq=23,
                 merge_distance=None,
                 threads=1):
        super().__init__()

        self._transposon_path = transposon_path
        self._index_path = bowtie_index_path
        self._bowtie_options = bowtie_options or {}

        self._min_length = min_length
        self._min_support = min_support
        self._min_mapq = min_mapq

        self._merge_distance = merge_distance
        self._threads = threads

    @classmethod
    def configure_args(cls, parser):
        cls._setup_base_args(parser, paired=True)

        parser.add_argument('--transposon', type=Path, required=True)
        parser.add_argument('--bowtie_index', type=Path, required=True)

        parser.add_argument('--min_length', type=int, default=15)
        parser.add_argument('--min_support', type=int, default=2)
        parser.add_argument('--min_mapq', type=int, default=23)
        parser.add_argument('--merge_distance', type=int, default=None)

        parser.add_argument('--local', default=False, action='store_true')
        parser.add_argument('--threads', default=1, type=int)

    @classmethod
    def _extract_args(cls, args):
        bowtie_options = {'--local': args.local, '--threads': args.threads}

        return dict(
            transposon_path=args.transposon,
            bowtie_index_path=args.bowtie_index,
            min_length=args.min_length,
            min_support=args.min_support,
            min_mapq=args.min_mapq,
            merge_distance=args.merge_distance,
            bowtie_options=bowtie_options,
            threads=args.threads)

    def run(self, read_path, output_dir, read2_path=None):
        if read2_path is None:
            raise ValueError('This pipeline requires paired-end data')

        logger = logging.getLogger()

        output_dir.mkdir(exist_ok=True, parents=True)

        # Extract genomic sequences.
        if logger is not None:
            logger.info('Extracting genomic sequences')
            logger.info('  %-18s: %s', 'Transposon',
                        shorten_path(self._transposon_path))
            logger.info('  %-18s: %s', 'Minimum length', self._min_length)

        # Trim reads and align to reference.
        genomic_path, genomic2_path = self._extract_genomic(
            read_path, read2_path, output_dir, logger=logger)
        alignment_path = self._align(
            genomic_path, genomic2_path, output_dir, logger=logger)

        # Extract insertions from bam file.
        bam_file = pysam.AlignmentFile(str(alignment_path))

        try:
            insertions = extract_insertions(
                iter(bam_file),
                func=_process_mates,
                paired=True,
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

    def _extract_genomic(self, read_path, read2_path, output_dir, logger):
        """Extracts genomic sequences from reads.

        Extracts genomic sequences by first trimming for mates for the
        transposon sequence (dropping reads without a match) and then
        trimming any Nextera transposase sequences from the remaining reads.
        Filtering for minimum length is performed in the nextera trimming step.
        """

        trimmed_tr_path, trimmed_tr2_path = self._trim_transposon(
            read_path, read2_path, output_dir, logger=logger)

        trimmed_nt_path, trimmed_nt2_path = self._trim_nextera(
            trimmed_tr_path, trimmed_tr2_path, output_dir, logger=logger)

        trimmed_tr_path.unlink()
        trimmed_tr2_path.unlink()

        return trimmed_nt_path, trimmed_nt2_path

    def _trim_transposon(self, read_path, read2_path, output_dir, logger):
        """Selects and trims mates with transposon sequence in second read."""

        cutadapt_opts = {'-G': 'file:' + str(self._transposon_path),
                         '--discard-untrimmed': True,
                         '--pair-filter=both': True}

        suffix = extract_suffix(read_path)
        trimmed_path = output_dir / ('genomic.R1' + suffix)
        trimmed2_path = output_dir / ('genomic.R2' + suffix)

        process = cutadapt(
            read_path=read_path,
            read2_path=read2_path,
            out_path=trimmed_path,
            out2_path=trimmed2_path,
            options=cutadapt_opts)

        if logger is not None:
            summary = cutadapt_summary(process.stdout, padding='   ')
            logger.info('Trimmed transposon sequence' + summary)

        return trimmed_path, trimmed2_path

    def _trim_nextera(self, read_path, read2_path, output_dir, logger):
        """Trims nextera sequences from mates and filters for min length."""

        cutadapt_opts = {
            '-a': 'CTGTCTCTTATA',
            '-A': 'CTGTCTCTTATA',
            '--minimum-length': self._min_length,
        }

        suffix = extract_suffix(read_path)
        trimmed_path = output_dir / ('trimmed_nextera.R1' + suffix)
        trimmed2_path = output_dir / ('trimmed_nextera.R2' + suffix)

        process = cutadapt(
            read_path=read_path,
            read2_path=read2_path,
            out_path=trimmed_path,
            out2_path=trimmed2_path,
            options=cutadapt_opts)

        if logger is not None:
            summary = cutadapt_summary(process.stdout, padding='    ')
            logger.info('Trimmed nextera sequences and '
                        'filtered for length' + summary)

        return trimmed_path, trimmed2_path

    def _align(self, read_path, read2_path, output_dir, logger):
        """Aligns mates to reference using bowtie2."""

        extra_opts = {'--threads': self._threads}
        options = toolz.merge(self._bowtie_options, extra_opts)

        # Align reads to genome.
        logger.info('Aligning to reference')
        logger.info('  %-18s: %s', 'Reference', shorten_path(self._index_path))
        logger.info('  %-18s: %s', 'Bowtie options',
                    flatten_arguments(options))

        alignment_path = output_dir / 'alignment.bam'

        bowtie2(
            read_paths=[read_path],
            read2_paths=[read2_path],
            index_path=self._index_path,
            output_path=alignment_path,
            options=options,
            verbose=True)

        return alignment_path


register_pipeline(name='nextera', pipeline=NexteraPipeline)


def _process_mates(mate1, mate2):
    ref = mate1.reference_name

    if mate1.is_reverse:
        transposon_pos = mate2.reference_start
        linker_pos = mate1.reference_end
        strand = 1
    else:
        transposon_pos = mate2.reference_end
        linker_pos = mate1.reference_start
        strand = -1

    return (ref, transposon_pos, strand), linker_pos
