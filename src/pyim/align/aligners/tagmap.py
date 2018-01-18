"""Module containing the nextera pipeline."""

import logging
from pathlib import Path

from pyim.model import InsertionSet
from pyim.util.path import extract_extension

from .base import Aligner

from ..base import (trim_transposon, trim_contaminants, align_reads,
                    extract_insertions)


class TagmapAligner(Aligner):
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
                 transposon_seq,
                 bowtie_index_path,
                 bowtie_options=None,
                 min_length=15,
                 min_support=2,
                 min_mapq=23,
                 merge_distance=None,
                 threads=1):
        super().__init__()

        self._transposon_seq = transposon_seq
        self._index_path = bowtie_index_path
        self._bowtie_options = bowtie_options or {}

        self._min_length = min_length
        self._min_support = min_support
        self._min_mapq = min_mapq

        self._merge_distance = merge_distance
        self._threads = threads

    def extract_insertions(self,
                           read_paths,
                           output_prefix=None,
                           verbose=False,
                           logger=None):
        """Extracts insertions from given reads."""

        output_prefix = Path(output_prefix)

        if logger is None:
            logger = logging.getLogger()

        if len(read_paths) != 2:
            raise ValueError('Pipeline only supports pair-end sequencing, '
                             'read_paths should contain 2 file paths')

        # Trim transposon sequence.
        if verbose:
            logger.info('- Trimming transposon sequences')

        read_ext = extract_extension(read_paths[0])
        trimmed_tr_paths = (
            Path(str(output_prefix) + '.trimmed_tr.R1' + read_ext),
            Path(str(output_prefix) + '.trimmed_tr.R2' + read_ext))

        trim_transposon(
            input_paths=read_paths,
            output_paths=trimmed_tr_paths,
            sequences={'-G': self._transposon_seq},
            threads=self._threads,
            verbose=verbose,
            logger=logger)

        # Remove nextera sequences.
        if verbose:
            logger.info('- Trimming nextera sequences')

        trimmed_nt_paths = (
            Path(str(output_prefix) + '.trimmed_nt.R1' + read_ext),
            Path(str(output_prefix) + '.trimmed_nt.R2' + read_ext))

        trim_contaminants(
            input_paths=trimmed_tr_paths,
            output_paths=trimmed_nt_paths,
            sequences={
                '-a': 'CTGTCTCTTATA',
                '-A': 'CTGTCTCTTATA',
            },
            min_length=self._min_length,
            threads=self._threads,
            verbose=verbose,
            logger=logger)

        # Align to reference.
        if verbose:
            logger.info('- Aligning to reference genome')

        alignment_path = Path(str(output_prefix) + '.bam')

        align_reads(
            read_paths=trimmed_nt_paths,
            index_path=self._index_path,
            output_path=alignment_path,
            threads=self._threads,
            verbose=verbose,
            logger=logger)

        # Extract insertions.
        if verbose:
            logger.info('- Extracting insertions')

        sample_name = output_prefix.name

        insertions = extract_insertions(
            alignment_path=alignment_path,
            position_func=self._position_for_mates,
            sample_func=lambda m1, m2: sample_name,
            min_mapq=self._min_mapq,
            paired=True,
            merge_dist=self._merge_distance,
            min_support=self._min_support)

        return InsertionSet.from_tuples(insertions)

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
