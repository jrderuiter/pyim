#! /usr/bin/env python

import os, argparse, logging, pandas
from os import path

from pyim.io import read_fasta
from pyim.ggplot import ggplot_save
from pyim.alignment.aligners import ExactReadAligner, TruncatedTargetAligner, ExonerateReadAligner, \
                                    ChainedReadAligner, ParallelTargetAligner, ParallelReadAligner
from pyim.alignment.filters import QueryEndFilter, MismatchFilter, BestScoreFilter
from pyim.vector.sb import SBVectorAligner, SBAlignmentStats
from pyim.genomic.bowtie2 import Bowtie2Aligner
from pyim.insertions import map_insertions, cluster_insertions

logging.basicConfig(format='%(asctime)s %(name)s \t %(message)s', level=logging.INFO)


def main(reads_file, reference, vector_file, barcode_file, barcode_mapping,
         contaminant_file, min_seq_length, output_dir, ncpus):
    # Create output dir if needed
    if not path.exists(output_dir):
        os.makedirs(output_dir)

    # Load reads
    logging.info('Reading fasta input')
    reads = list(read_fasta(reads_file))

    # Remove contaminated reads
    if contaminant_file is not None:
        logging.info('Removing contaminated reads')
        num_reads = len(reads)
        for contaminant in read_fasta(contaminant_file):
           reads = [r for r in reads if contaminant.seq not in r.seq]
        num_dropped =  num_reads - len(reads)
        logging.info('  Dropped %d (%3.2f%%) contaminated reads' % (num_dropped, (num_dropped / float(num_reads))*100))

    # Do vector alignment (align SB, T7 and barcodes)
    logging.info('Aligning vectors to reads')
    vec_aln = _align_vectors(reads, vector_file, barcode_file, ncpus)
    _vec_stats(reads, vec_aln, output_dir)

    logging.info('Filtering short genomic sequences')
    long_seqs = [s for s in vec_aln.genomic_sequences if len(s.seq) >= min_seq_length]

    # Do genomic alignment
    logging.info('Aligning genomic sequences')
    gen_aligner = Bowtie2Aligner(num_cores=ncpus)
    hits, uniq, dupl, unaligned = gen_aligner.align(long_seqs, reference)

    from pyim.io import write_fasta
    write_fasta('unaligned.fa', unaligned)

    # Determine and write insertions
    logging.info('Mapping and clustering putative insertions')
    insertions_raw = map_insertions(hits, vec_aln.barcode_alignments)
    insertions = cluster_insertions(insertions_raw, dist_t=10)

    # Add sample names
    if barcode_mapping is not None:
        logging.info('Looking up sample names from barcode mapping')
        bc_map = pandas.read_csv(barcode_mapping)
        bc_dict = dict(zip(bc_map['Barcode'], bc_map['Sample']))

        insertions_raw['sample'] = insertions_raw['barcode'].map(bc_dict)
        insertions['sample'] = insertions['barcode'].map(bc_dict)

    # Write output
    logging.info('Writing output')
    insertions.to_csv('insertions.csv', sep='\t')


def _align_vectors(reads, vector_file, barcode_file, ncpus):
    trunc1_aligner = TruncatedTargetAligner(ExactReadAligner(), 1)
    trunc3_aligner = TruncatedTargetAligner(ExactReadAligner(), 3)
    ex_aligner = ExonerateReadAligner(min_score=60, filters=[QueryEndFilter(verbose=True)])

    exhaust_t7_filters = [QueryEndFilter(verbose=True),
                          MismatchFilter(2, include_start=True, verbose=True),
                          BestScoreFilter(verbose=True)]
    exhaust_aligner = ExonerateReadAligner(min_score=40, exhaustive=True, bestn=10)
    exhaust_aligner = ParallelReadAligner(exhaust_aligner, ncpus, filters=exhaust_t7_filters)

    t7_aligner = ChainedReadAligner([ExactReadAligner(), trunc1_aligner, trunc3_aligner, ex_aligner, exhaust_aligner])
    bc_aligner = ParallelTargetAligner(ExactReadAligner(), ncpus)
    sb_aligner = ChainedReadAligner([ExactReadAligner(), ExonerateReadAligner(min_score=60)])

    vec_aligner = SBVectorAligner()
    aln_result = vec_aligner.align(reads, vector_file, barcode_file,
                                   t7_aligner=t7_aligner, bc_aligner=bc_aligner, sb_aligner=sb_aligner)

    return aln_result

def _vec_stats(reads, aln_result, output_dir):
    stats = SBAlignmentStats(reads, aln_result)
    ggplot_save(stats.plot_mapped(),                  path.join(output_dir, 'vectors_mapped.pdf'))
    ggplot_save(stats.plot_unmapped_type(kind='bar'), path.join(output_dir, 'vectors_unmapped_type.pdf'))
    ggplot_save(stats.plot_alignment_types(),         path.join(output_dir, 'vectors_alignment_type.pdf'))


def _parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument('--reads', required=True)
    parser.add_argument('--reference', required=True)

    parser.add_argument('--output_dir', default='.')

    parser.add_argument('--vector_file', required=True)
    parser.add_argument('--barcode_file', required=True)
    parser.add_argument('--barcode_mapping', default=None)
    parser.add_argument('--contaminant_file', default=None)

    parser.add_argument('--min_seq_len', default=15, type=int)
    parser.add_argument('--ncpus', default=1, type=int)

    return parser.parse_args()



if __name__ == '__main__':
    args = _parse_args()
    main(args.reads, args.reference, args.vector_file, args.barcode_file,
         args.barcode_mapping, args.contaminant_file, args.min_seq_len, args.output_dir, args.ncpus)

