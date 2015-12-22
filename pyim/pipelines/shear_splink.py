import os
import operator
import logging
from enum import Enum
from os import path

import pysam
import pandas as pd
import skbio
import toolz
import tqdm
from toolz.curried import (filter as curried_filter,
                           map as curried_map)

from pyim.alignment.bowtie2 import align as bowtie_align
from pyim.alignment.vector import (align_exact, align_multiple,
                                   align_with_reverse, reverse_alignment)
from pyim.util import count_fasta_entries

from ._model import ExtractResult
from ._helpers.pipeline import (print_stats, build_barcode_map,
                                write_genomic_sequences)
from ._helpers.grouping import (chain_groupby, groupby_barcode,
                                groupby_reference_position)
from ._helpers.clustering import merge_within_distance


logging.basicConfig(
        format='%(asctime)-15s %(message)s',
        datefmt='[%Y-%m-%d %H:%M:%S]',
        level=logging.INFO)


# --- Pipeline register hook + main --- #

def register(subparsers, name='shear_splink'):
    parser = subparsers.add_parser(name, help=name + ' help')

    # Required arguments.
    parser.add_argument('input')
    parser.add_argument('output_dir')
    parser.add_argument('bowtie_index')
    parser.add_argument('transposon')
    parser.add_argument('barcodes')
    parser.add_argument('linker')

    # Optional arguments.
    parser.add_argument('--contaminants', default=None)
    parser.add_argument('--sample_map', default=None)
    parser.add_argument('--min_genomic_length', type=int, default=15)
    parser.add_argument('--min_depth', type=int, default=2)
    parser.add_argument('--min_mapq', type=int, default=37)

    # Set main for dispatch.
    parser.set_defaults(main=main)

    return parser


def main(args):
    # Read transposon, linker and barcode sequences.
    transposon = skbio.io.read(args.transposon, format='fasta', into=skbio.DNA)
    linker = skbio.io.read(args.linker, format='fasta', into=skbio.DNA)

    barcodes = list(skbio.io.read(args.barcodes, format='fasta',
                                  constructor=skbio.DNA))

    if args.contaminants is not None:
        contaminants = list(skbio.io.read(args.contaminants, format='fasta',
                                          constructor=skbio.DNA))
    else:
        contaminants = None

    # Read barcode --> sample map if given.
    if args.sample_map is not None:
        sample_map = pd.read_csv(args.sample_map, sep='\t')
        sample_map = dict(zip(sample_map['barcode'],
                              sample_map['sample']))
    else:
        sample_map = None

    # Create output_dir if it does not exist.
    if not path.exists(args.output_dir):
        os.mkdir(args.output_dir)

    # Run pipeline!
    insertions = shear_splink(
        args.input, transposon, linker, barcodes,
        args.bowtie_index, args.output_dir,
        contaminants=contaminants, sample_map=sample_map,
        min_genomic_length=args.min_genomic_length)

    # Write insertion output.
    insertions.to_csv(path.join(args.output_dir, 'insertions.txt'),
                      sep='\t', index=False)


# --- Overall pipeline --- #

def shear_splink(read_path, transposon, linker, barcodes,
                 bowtie_index, output_dir, contaminants=None,
                 sample_map=None, min_genomic_length=15):

    logger = logging.getLogger()

    # Determine paths for intermediates/outputs.
    genomic_path = path.join(output_dir, 'genomic.fna')
    barcode_path = path.join(output_dir, 'genomic.barcodes.txt')
    alignment_base = path.join(output_dir, 'alignment')

    # Log progress with progressbar.
    logger.info('Extracting genomic sequences')

    reads = skbio.read(read_path, format='fasta')
    reads = tqdm.tqdm(reads, total=count_fasta_entries(read_path),
                      leave=False, ncols=80)

    # Extract genomic sequences and barcodes
    _, barcode_frame = extract_genomic(
        reads, transposon=transposon, barcodes=barcodes, linker=linker,
        output_path=genomic_path, contaminants=contaminants,
        min_length=min_genomic_length, logger=logger)

    barcode_frame.to_csv(barcode_path, sep='\t', index=False)

    # Align to reference with Bowtie2.
    logger.info('Aligning to reference genome')

    aln_path = bowtie_align(genomic_path, bowtie_index, alignment_base,
                            bam_output=True, options={'-f': True},
                            log=alignment_base + '.log')

    # Identify insertions from alignment.
    logger.info('Identifying insertions')

    barcode_map = dict(zip(barcode_frame['read_id'],
                           barcode_frame['barcode']))
    insertions = identify_insertions(aln_path, barcode_map=barcode_map)

    # Cluster and merge close insertions
    logger.info('Merging close insertions')

    agg_funcs = {'depth': 'sum', 'depth_unique': 'sum'}
    insertions = merge_within_distance(
        insertions, max_dist=2000, agg_funcs=agg_funcs)

    # Assign ids to insertions.
    insertions['id'] = ['INS_{}'.format(i)
                        for i in range(1, len(insertions) + 1)]

    # Map barcodes to samples.
    if sample_map is not None:
        insertions['sample'] = insertions['barcode'].map(sample_map)

    return insertions


# --- Genomic sequence extraction --- #

class ShearSplinkStatus(Enum):
    contaminant = 1
    no_transposon = 2
    no_linker = 3
    no_barcode = 4
    multiple_barcodes = 5
    too_short = 6
    proper_read = 7


def extract_genomic(reads, transposon, barcodes, linker,
                    output_path, contaminants=None, min_length=15,
                    io_kwargs=None, logger=None):
    io_kwargs = io_kwargs or {}

    # Extract and write genomic sequences.
    barcode_map = toolz.pipe(
        reads,
        _extract_from_reads(transposon=transposon,
                            barcodes=barcodes,
                            linker=linker,
                            contaminants=contaminants),
        curried_map(_check_minimum_length(min_length=15)),
        print_stats(logger=logger),
        curried_filter(_proper_filter),
        write_genomic_sequences(file_path=output_path,
                                format='fasta', **io_kwargs),
        build_barcode_map)

    # Build frame mapping reads to barcodes.
    barcode_frame = pd.DataFrame.from_records(
        iter(barcode_map.items()), columns=['read_id', 'barcode'])

    return output_path, barcode_frame


@toolz.curry
def _extract_from_reads(
        reads, transposon, barcodes, linker, contaminants=None,
        transposon_func=None, barcode_func=None, linker_func=None):

    # Specify defaults for not provided aligners.
    if transposon_func is None:
        transposon_func = align_with_reverse(align_func=align_exact)

    if barcode_func is None:
        barcode_func = align_multiple(align_func=align_exact)

    if linker_func is None:
        linker_func = align_exact

    # Setup contaminant aligner if sequences are provided.
    if contaminants is not None:
        contaminant_func = align_multiple(queries=contaminants,
                                          align_func=align_exact,
                                          return_first=True)
    else:
        contaminant_func = None

    # Prime aligners with their respective sequences.
    transposon_func = transposon_func(query=transposon)
    barcode_func = barcode_func(queries=barcodes)
    linker_func = linker_func(query=linker)

    # Extract and return results.
    extract_func = toolz.curry(
         _extract_from_read,
         transposon_func=transposon_func,
         barcode_func=barcode_func,
         linker_func=linker_func,
         contaminant_func=contaminant_func)

    for result in map(extract_func, reads):
        yield result


def _extract_from_read(read, transposon_func, barcode_func,
                       linker_func, contaminant_func=None):
    """ Extracts the genomic sequence and barcode from the passed
        read. Reads containing contaminants are dropped. Reads are
        expected to look as follows:

            [barcode][transposon][genomic-sequence][linker]

        Each of these sequences is recognized by their corresponding
        alignment function. The barcode alignment identifies the
        barcode (and thus the sample) of the read, whilst the transposon
        and linker alignments are used to delineate the genomic sequence.

        The function returns an ExactResult tuple that contains the
        genomic sequence, barcode and a status flag. If any errors
        occur during the extraction, the genomic sequence and barcode
        values are None and the status flag indicates the underlying reason.
    """

    # Drop read if it contains a contaminant.
    if contaminant_func is not None and contaminant_func(read) is not None:
        return ExtractResult(None, None, ShearSplinkStatus.contaminant)

    # Identify location of the transposon.
    transposon_aln = transposon_func(read)
    if transposon_aln is None:
        return ExtractResult(None, None, ShearSplinkStatus.no_transposon)

    # If transposon is on the reverse strand, flip the read and the
    # alignment to bring everything into the same (fwd) orientation.
    if transposon_aln.strand == -1:
        read = read.reverse_complement()
        transposon_aln = reverse_alignment(transposon_aln)

    # Identify barcode of the read.
    try:
        barcode_aln = barcode_func(read)
        if barcode_aln is None:
            return ExtractResult(None, None, ShearSplinkStatus.no_barcode)
    except ValueError:
        return ExtractResult(None, None, ShearSplinkStatus.multiple_barcodes)

    barcode = barcode_aln.query_id

    # Identify location of linker.
    linker_aln = linker_func(read)
    if linker_aln is None:
        return ExtractResult(None, None, ShearSplinkStatus.no_linker)

    # Extract genomic sequence using previous alignments.
    genomic = read[transposon_aln.target_end:linker_aln.target_start]

    return ExtractResult(genomic, barcode, ShearSplinkStatus.proper_read)


@toolz.curry
def _check_minimum_length(result, min_length):
    """Flags proper reads if shorter than min_length."""
    if (result.status == ShearSplinkStatus.proper_read and
            len(result.genomic_sequence) < min_length):
        result.status = ShearSplinkStatus.too_short
    return result


def _proper_filter(result):
    """Filters extraction results for proper reads."""
    return result.status == ShearSplinkStatus.proper_read


# --- Insertion identification --- #

def identify_insertions(alignment_path, barcode_map):
    # Get alignments from bowtie.
    bam = pysam.AlignmentFile(alignment_path)
    alignments = bam.fetch(multiple_iterators=True)

    # Group alignments by barcode and position.
    aln_groups = chain_groupby(
        alignments,
        [groupby_reference_position(alignment_file=bam),
         groupby_barcode(barcode_map=barcode_map)])

    # Convert groups into insertion frame.
    insertions = pd.DataFrame.from_records(
        (_alignments_to_insertion(info, alns)
         for info, alns in aln_groups) ,
        columns=['id', 'chrom', 'position', 'strand',
                 'barcode', 'depth', 'depth_unique'])

    return insertions


def _alignments_to_insertion(info, alignments, id_=None):
    # Extract group info.
    ref, pos, strand, bc = info

    # Get positions of the non-transposon ends of the alignment.
    end_field = 'reference_end' if strand == 1 else 'reference_start'
    end_positions = map(operator.attrgetter(end_field), alignments)

    # Calculate overall depth and unique end depth.
    depth = len(alignments)
    depth_unique = len(set(end_positions))

    return id_, ref, pos, strand, bc, depth, depth_unique
