import logging
from pathlib import Path

import pysam

from pyim.external.bowtie2 import bowtie2
from pyim.external.cutadapt import cutadapt
from pyim.external.util import flatten_arguments
from pyim.model import Insertion
from pyim.util.path import shorten_path

from ..common.insertions import extract_insertions
from .base import Pipeline


class NexteraPipeline(Pipeline):
    """Nextera-based transposon pipeline."""

    def __init__(self,
                 transposon_path,
                 read1_adapter_path,
                 read2_adapter_path,
                 min_length=20,
                 logger=None):
        super().__init__(logger=logger)

        self._read1_adapter_path = read1_adapter_path
        self._read2_adapter_path = read2_adapter_path
        self._transposon_path = transposon_path
        self._min_length = min_length

    def _extract_args(cls, args):
        raise NotImplementedError()

    def run(self, read_path, output_dir, read2_path=None):
        logger = logging.getLogger()

        trimmed_path, trimmed2_path = self._trim_adapters(
            read_path, read2_path, output_dir, logger=logger)
        alignment_path = self._align(
            trimmed_path, trimmed2_path, output_dir, logger=logger)

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

        # genomic_paths = self._extract_genomic(trimmed_paths, work_dir)
        # alignment_path = self._align(genomic_paths)
        # yield from self._extract_insertions(alignment_path)

        # def _trim_nextera_adapters(self, read_paths, work_dir):
        #     output_paths = path.build_paths(
        #         read_paths, dir_=work_dir / '_trim_nextera')

        #     cutadapt_opts = {'--minimum-length': self._min_length,
        #                      '-g': 'file:' + str(self._read1_adapter_path),
        #                      '-G': 'file:' + str(self._read2_adapter_path)}

        #     cutadapt(
        #         read_paths[0],
        #         output_paths[0],
        #         cutadapt_opts,
        #         read2_path=read_paths[1],
        #         out2_path=output_paths[1])

        #     return output_paths

        # def _extract_genomic(self, read_paths, work_dir):
        #     output_paths = path.build_paths(
        #         read_paths, dir_=work_dir / '_trim_sequence')

        #     cutadapt_opts = {'--minimum-length': self._min_length,
        #                      '--discard-untrimmed': True,
        #                      '--pair-filter=both': True,
        #                      '-G': 'file:' + str(self._transposon_path)}

        #     cutadapt(
        #         read_paths[0],
        #         output_paths[0],
        #         cutadapt_opts,
        #         read2_path=read_paths[1],
        #         out2_path=output_paths[1])

        #     return output_paths

    def _trim_adapters(self, read_path, read2_path, output_dir, logger):
        raise NotImplementedError()

    def _align(self, read_path, read2_path, output_dir, logger):

        # Align reads to genome.
        logger.info('Aligning to reference')
        logger.info('  %-18s: %s', 'Reference', shorten_path(self._index_path))
        logger.info('  %-18s: %s', 'Bowtie options',
                    flatten_arguments(self._bowtie_options))

        alignment_path = output_dir / 'alignment.bam'

        bowtie2(
            read_paths=[read_path],
            read2_paths=[read2_path],
            index_path=self._index_path,
            output_path=alignment_path,
            options=self._bowtie_options,
            verbose=True)

        return alignment_path


def _trim_adapters(read_paths,
                   output_paths,
                   adapter1_path=None,
                   adapter2_path=None,
                   discard=False,
                   min_length=None):
    cutadapt_opts = {}

    if adapter1_path is not None:
        cutadapt_opts['-g'] = 'file:' + str(adapter1_path)

    if adapter2_path is not None:
        cutadapt_opts['-G'] = 'file:' + str(adapter2_path)

    if discard:
        cutadapt_opts['--discard-trimmed'] = True

    if min_length is not None:
        cutadapt_opts['--minimum-length'] = min_length

    cutadapt(
        read_paths[0],
        output_paths[0],
        cutadapt_opts,
        in2_path=read_paths[1],
        out2_path=output_paths[1])


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
