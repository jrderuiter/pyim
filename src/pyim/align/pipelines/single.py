import itertools
import logging
import os
from pathlib import Path

from cutadapt import seqio
import pandas as pd

from pyim.external import bowtie2
from pyim.external.util import flatten_options
from pyim.util.path import build_path

from ..common import genomic as cm_gen, insertions as cm_ins
from .base import Pipeline, register_pipeline


class SinglePipeline(Pipeline):
    def __init__(self,
                 transposon_path,
                 bowtie_index_path,
                 linker_path=None,
                 contaminant_path=None,
                 min_length=15,
                 min_support=2,
                 min_mapq=23,
                 merge_distance=0,
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
        super().configure_args(parser)

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
    def extract_args(cls, args):
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

    def run(self, reads_path, work_dir):
        logger = logging.getLogger()

        # Extract genomic sequences and align to reference.
        alignment_path = self._extract_and_align(reads_path, work_dir, logger)

        # Extract alignment groups (grouped by position) from bam file.
        logger.info('Summarizing alignments')
        logger.info('  %-18s: %s', 'Minimum mapq', self._min_mapq)

        alignments = cm_ins.fetch_alignments(
            alignment_path, min_mapq=self._min_mapq)
        aln_summary = cm_ins.summarize_alignments(alignments)

        # Convert groups to insertions and return.
        logger.info('Converting to insertions')
        logger.info('  %-18s: %d', 'Minimum support', self._min_support)
        logger.info('  %-18s: %d', 'Merge distance', self._merge_distance)

        yield from cm_ins.convert_groups_to_insertions(
            aln_summary,
            min_support=self._min_support,
            merge_distance=self._merge_distance)

    def _extract_and_align(self, reads_path, work_dir, logger):
        # Extract genomic sequences.
        logger.info('Extracting genomic sequences')
        logger.info('  %-18s: %s', 'Transposon',
                    shorten_path(self._transposon_path))
        logger.info('  %-18s: %s', 'Linker', shorten_path(self._linker_path))
        logger.info('  %-18s: %s', 'Contaminants',
                    shorten_path(self._contaminant_path))
        logger.info('  %-18s: %s', 'Minimum length', self._min_length)

        genomic_path = build_path(reads_path, dir_=work_dir, suffix='.genomic')
        genomic_path.parent.mkdir(exist_ok=True, parents=True)

        cm_gen.extract_genomic(
            reads_path,
            genomic_path,
            transposon_path=self._transposon_path,
            linker_path=self._linker_path,
            contaminant_path=self._contaminant_path,
            min_length=self._min_length,
            min_overlaps=self._min_overlaps,
            error_rates=self._error_rates)

        # Align reads to genome.
        logger.info('Aligning to reference')
        logger.info('  %-18s: %s', 'Reference', shorten_path(self._index_path))
        logger.info('  %-18s: %s', 'Bowtie options',
                    flatten_options(self._bowtie_options))

        alignment_path = build_path(reads_path, dir_=work_dir, ext='.bam')
        alignment_path.parent.mkdir(exist_ok=True, parents=True)

        bowtie2.bowtie2(
            [genomic_path],
            self._index_path,
            alignment_path,
            options=self._bowtie_options,
            verbose=True)

        return alignment_path


register_pipeline(name='single', pipeline=SinglePipeline)


class SingleMultiplexedPipeline(SinglePipeline):
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
    def extract_args(cls, args):
        arg_dict = super().extract_args(args)

        if args.barcode_mapping is not None:
            map_df = pd.read_csv(args.barcode_mapping, sep='\t')
            arg_dict['barcode_mapping'] = dict(
                zip(map_df['barcode'], map_df['sample']))
        else:
            arg_dict['barcode_mapping'] = None

        arg_dict['barcode_path'] = args.barcodes

        return arg_dict

    def run(self, reads_path, work_dir):
        logger = logging.getLogger()

        # Extract genomic sequences and align to reference.
        alignment_path = self._extract_and_align(reads_path, work_dir, logger)

        # Map reads to specific barcodes/samples.
        logger.info('Extracting barcode/sample mapping')
        logger.info('  %-18s: %s', 'Barcodes',
                    shorten_path(self._barcode_path))
        read_map = self._get_barcode_mapping(reads_path)

        # Extract alignment groups (grouped by position) from bam file.
        logger.info('Summarizing alignments')
        logger.info('  %-18s: %s', 'Minimum mapq', self._min_mapq)

        alignments = cm_ins.fetch_alignments(
            alignment_path, min_mapq=self._min_mapq)

        aln_summaries = cm_ins.summarize_alignments_by_group(
            alignments,
            group_func=lambda aln: read_map.get(aln.query_name, None))

        # Convert groups from each sample into insertions,
        # adding sample name and sample prefix to the ID.
        logger.info('Converting to insertions')
        logger.info('  %-18s: %d', 'Minimum support', self._min_support)
        logger.info('  %-18s: %d', 'Merge distance', self._merge_distance)

        insertion_grps = (
            cm_ins.convert_summary_to_insertions(
                aln_summ,
                min_support=self._min_support,
                merge_distance=self._merge_distance,
                sample=barcode,
                id_fmt=barcode + '.INS_{}')
            for barcode, aln_summ in aln_summaries.items()) # yapf: disable

        # Return concatenated list of insertions.
        yield from itertools.chain.from_iterable(insertion_grps)

    def _get_barcode_mapping(self, reads_path):
        # Read barcode sequences.
        with seqio.open(str(self._barcode_path)) as barcode_file:
            barcodes = list(barcode_file)

        # Extract read --> barcode mapping.
        with seqio.open(str(reads_path)) as reads:
            return cm_ins.extract_barcode_mapping(reads, barcodes,
                                                  self._barcode_mapping)


register_pipeline(
    name='single-multiplexed', pipeline=SingleMultiplexedPipeline)


def shorten_path(file_name, limit=40):
    f = os.path.split(str(file_name))[1]
    return "%s~%s" % (f[:3], f[-(limit - 3):]) if len(f) > limit else f
