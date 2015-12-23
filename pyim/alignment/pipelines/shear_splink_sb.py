import os
from os import path

import pandas as pd
import skbio
from .shear_splink import shear_splink

from pyim.alignment import vector as vec
from pyim.util import count_fasta_entries


# --- Pipeline register hook + main --- #

def register(subparsers, name='shear_splink_sb'):
    parser = subparsers.add_parser(name, help=name + ' help')

    # Required arguments.
    parser.add_argument('input')
    parser.add_argument('output_dir')
    parser.add_argument('--bowtie_index', required=True)
    parser.add_argument('--transposon', required=True)
    parser.add_argument('--barcodes', required=True)
    parser.add_argument('--linker', required=True)

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
    # Prepare reads, counting total for progress bar.
    reads = skbio.read(args.input, format='fasta', constructor=skbio.DNA)
    total_reads = count_fasta_entries(args.input)

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
        os.makedirs(args.output_dir, exist_ok=True)

    # Setup custom aligners.
    transposon_aligner = vec.align_chained(
        align_funcs=[vec.compose(vec.align_exact, try_reverse=True),
                     vec.compose(vec.align_ssw, try_reverse=True,
                                 filters=[vec.filter_score(min_score=90)])])

    linker_ssw_filters = [
        vec.filter_score(min_score=90),
        vec.filter_and(filters=[
            vec.filter_end_match(),
            vec.filter_coverage(min_coverage=0.5, min_identity=0.9)])]

    linker_aligner = vec.align_chained(
        align_funcs=[vec.compose(vec.align_exact, try_reverse=True),
                     vec.compose(vec.align_ssw, try_reverse=True,
                                 filters=linker_ssw_filters)])

    extract_kws = {'linker_func': linker_aligner,
                   'transposon_func': transposon_aligner}

    # Run pipeline!
    insertions = shear_splink(
        reads, transposon, linker, barcodes,
        args.bowtie_index, args.output_dir,
        contaminants=contaminants, sample_map=sample_map,
        min_genomic_length=args.min_genomic_length,
        extract_kws=extract_kws, total_reads=total_reads)

    # Write insertion output.
    insertions.to_csv(path.join(args.output_dir, 'insertions.txt'),
                      sep='\t', index=False)
