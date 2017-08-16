"""Module containing the nextera pipeline."""

import logging
from pathlib import Path

import pysam
import toolz

from pyim.external.bowtie2 import bowtie2
from pyim.external.cutadapt import cutadapt, cutadapt_summary
from pyim.external.util import flatten_arguments
from pyim.model import Insertion
from pyim.util.path import WorkDirectory, shorten_path, extract_suffix

from .base import Aligner, PairedEndCommand
from ..util import AlignmentSummary


class NexteraAligner(Aligner):
    """Nextera-based transposon pipeline.

    Analyzes paired-end sequence data that was prepared using a Nextera-based
    protocol. Sequence reads are expected to have the following structure::

        Mate 1:
            [Genomic]

        Mate 2:
            [Transposon][Genomic]

    Here, ``transposon`` refers to the flanking part of the transposon sequence
    and ``genomic`` refers to the genomic DNA located between the transposon
    sequence and the used adapt sequence. Note that the adapter itself is not
    sequenced and therefore not part of the reads. However, the end of Mate 1
    is considered to terminate at the adapter and as such represents the
    breakpoint between the genomic DNA and the adapter.

    The pipeline essentially performs the following steps:

        - Mates are trimmed to remove the transposon sequence, dropping any
          reads not containing the transposon.
        - The remaining mates are trimmed to remove any sequences from
          the Nextera transposase.
        - The trimmed mates are aligned to the reference genome.
        - The resulting alignment is used to identify insertions.

    Parameters
    ----------
    transposon_path : Path
        Path to the (flanking) transposon sequence (fasta).
    bowtie_index_path : Path
        Path to the bowtie index.
    bowtie_options : Dict[str, Any]
        Dictionary of extra options for Bowtie.
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
    threads : int
        The number of threads to use for the alignment.

    """

    def __init__(self,
                 transposon_path,
                 bowtie_index_path,
                 bowtie_options=None,
                 min_length=15,
                 min_support=2,
                 min_mapq=23,
                 merge_distance=None,
                 threads=1,
                 sample_name=None,
                 logger=None):
        super().__init__()

        self._transposon_path = transposon_path
        self._index_path = bowtie_index_path
        self._bowtie_options = bowtie_options or {}

        self._min_length = min_length
        self._min_support = min_support
        self._min_mapq = min_mapq

        self._merge_distance = merge_distance
        self._threads = threads

        self._sample_name = sample_name

        self._logger = logger or logging.getLogger()

    def trim(self, read_paths, output_paths, work_dir=None):
        """Trims reads to remove transposon/nextera sequences."""

        # if logger is not None:
        #     logger.info('Extracting genomic sequences')
        #     logger.info('  %-18s: %s', 'Transposon',
        #                 shorten_path(self._transposon_path))
        #     logger.info('  %-18s: %s', 'Minimum length', self._min_length)

        self._check_read_paths(read_paths)
        suffix = extract_suffix(read_paths[0])

        with WorkDirectory(work_dir, keep=work_dir is not None) as work_dir:
            # Select reads with transposon and trim sequence.
            trimmed_tr_paths = (work_dir / ('trimmed.transposon.R1' + suffix),
                                work_dir / ('trimmed.transposon.R2' + suffix))
            self._trim_transposon(read_paths, trimmed_tr_paths)

            # Trim nextera sequences if present.
            trimmed_nt_paths = (work_dir / ('trimmed.nextera.R1' + suffix),
                                work_dir / ('trimmed.nextera.R2' + suffix))
            self._trim_nextera(trimmed_tr_paths, trimmed_nt_paths)

            # Move outputs into position.
            for file_path, output_path in zip(trimmed_nt_paths, output_paths):
                file_path.rename(output_path)

    def _check_read_paths(self, read_paths):
        """Checks read paths input for validity."""

        if len(read_paths) != 2:
            raise ValueError(self.__class__.__name__ +
                             ' only supports paired-end data')

    def _trim_transposon(self, read_paths, output_paths):
        """Selects and trims mates with transposon sequence in second read."""

        cutadapt_opts = {
            '-G': 'file:' + str(self._transposon_path),
            '--discard-untrimmed': True,
            '--pair-filter=both': True
        }

        process = cutadapt(
            read_path=read_paths[0],
            read2_path=read_paths[1],
            out_path=output_paths[0],
            out2_path=output_paths[1],
            options=cutadapt_opts)

        summary = cutadapt_summary(process.stdout, padding='   ')
        self._logger.info('Trimmed transposon sequence' + summary)

    def _trim_nextera(self, read_paths, output_paths):
        """Trims nextera sequences from mates and filters for min length."""

        cutadapt_opts = {
            '-a': 'CTGTCTCTTATA',
            '-A': 'CTGTCTCTTATA',
            '--minimum-length': self._min_length,
        }

        process = cutadapt(
            read_path=read_paths[0],
            read2_path=read_paths[1],
            out_path=output_paths[0],
            out2_path=output_paths[1],
            options=cutadapt_opts)

        summary = cutadapt_summary(process.stdout, padding='    ')
        self._logger.info('Trimmed nextera sequences and '
                          'filtered for length' + summary)

    def align(self, read_paths, output_path):
        """Aligns mates to reference using bowtie2."""

        self._check_read_paths(read_paths)

        extra_opts = {'--threads': self._threads}
        options = toolz.merge(self._bowtie_options, extra_opts)

        # Align reads to genome.
        # logger.info('Aligning to reference')
        # logger.info('  %-18s: %s', 'Reference', shorten_path(self._index_path))
        # logger.info('  %-18s: %s', 'Bowtie options',
        #             flatten_arguments(options))

        bowtie2(
            read_paths=[read_paths[0]],
            read2_paths=[read_paths[1]],
            index_path=self._index_path,
            output_path=output_path,
            options=options,
            verbose=True)

    def extract(self, bam_path):
        """Extract insertions from alignment."""

        bam_file = pysam.AlignmentFile(str(bam_path))

        try:
            summary = AlignmentSummary.from_alignments(
                iter(bam_file),
                position_func=self._position_for_mates,
                sample_func=lambda m1, m2: self._sample_name,
                min_mapq=self._min_mapq,
                paired=True)
        finally:
            bam_file.close()

        if self._merge_distance is not None:
            summary = summary.merge_within_distance(self._merge_distance)

        insertions = summary.to_insertions(min_support=self._min_support)

        yield from insertions

    @staticmethod
    def _position_for_mates(mate1, mate2):
        """Returns transposon/linker positions for given mates."""

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

    def run(self, read_paths, work_dir=None):
        """Runs aligner on given read files."""

        self._check_read_paths(read_paths)
        suffix = extract_suffix(read_paths[0])

        with WorkDirectory(work_dir, keep=work_dir is not None) as work_dir:
            # Trim reads and align to reference.
            trimmed_paths = (work_dir / ('trimmed.R1' + suffix),
                             work_dir / ('trimmed.R2' + suffix))
            self.trim(read_paths, trimmed_paths, work_dir=work_dir)

            alignment_path = work_dir / 'alignment.bam'
            self.align(trimmed_paths, output_path=alignment_path)

            # Extract insertions.
            insertions = list(self.extract(alignment_path))

        return insertions


class NexteraCommand(PairedEndCommand):
    """Command for the Nextera aligner."""

    name = 'nextera'

    def configure(self, parser):
        super().configure(parser)

        parser.add_argument('--transposon', type=Path, required=True)
        parser.add_argument('--bowtie_index', type=Path, required=True)

        parser.add_argument('--sample_name', required=True)

        parser.add_argument('--min_length', type=int, default=15)
        parser.add_argument('--min_support', type=int, default=2)
        parser.add_argument('--min_mapq', type=int, default=23)
        parser.add_argument('--merge_distance', type=int, default=None)

        parser.add_argument('--local', default=False, action='store_true')

        parser.add_argument('--work_dir', default=None)
        parser.add_argument('--threads', default=1, type=int)

        return parser

    def run(self, args):
        bowtie_options = {'--local': args.local, '--threads': args.threads}

        aligner = NexteraAligner(
            transposon_path=args.transposon,
            bowtie_index_path=args.bowtie_index,
            min_length=args.min_length,
            min_support=args.min_support,
            min_mapq=args.min_mapq,
            merge_distance=args.merge_distance,
            bowtie_options=bowtie_options,
            sample_name=args.sample_name,
            threads=args.threads)

        insertions = aligner.run(args.reads, work_dir=args.work_dir)

        args.output.parent.mkdir(exist_ok=True, parents=True)
        Insertion.to_csv(args.output, insertions, sep='\t', index=False)
