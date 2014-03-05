
import logging

from pyim.io import read_fasta
from pyim.alignment.genomic.bowtie2 import Bowtie2Aligner
from pyim.insertions.clustering import cluster_insertions
from pyim.insertions.mapping import map_insertions_to_samples

from pyim.sb.alignment import align
from pyim.sb.insertions import hits_to_insertions


logging.basicConfig(format='%(asctime)s %(name)s \t %(message)s', level=logging.INFO)


def sb_pipeline(reads_file, reference, vector_file, barcode_file, barcode_sample_file,
                contaminant_file=None, threads=1, min_sequence_len=15, min_ulp=2,
                sb_aligner=None, t7_aligner=None, bc_aligner=None):

    logging.info('Reading input reads')
    reads = list(read_fasta(reads_file))

    logging.info('Filtering contaminants')
    if contaminant_file is not None:
        reads = filter_contaminated_reads(reads, contaminant_file)

    logging.info('Aligning vectors and extracting genomic sequences')
    genomic_seqs, barcode_mapping = align(reads, vector_file, barcode_file, sb_aligner=sb_aligner,
                                          t7_aligner=t7_aligner, bc_aligner=bc_aligner)

    logging.info('Filtering short reads')
    genomic_seqs = [g for g in genomic_seqs if len(g.seq) >= min_sequence_len]

    logging.info('Aligning genomic sequences to reference genome')
    genomic_aligner = Bowtie2Aligner(num_cores=threads)
    hits, _, _, _ = genomic_aligner.align(genomic_seqs, reference)

    # Derive insertions from hits and cluster close insertions in the
    # same sample to avoid double hits due to sequencing errors
    logging.info('Generating insertions from alignments')
    insertions = hits_to_insertions(hits, barcode_mapping)
    insertions = cluster_insertions(insertions, max_dist=10)

    logging.info('Restricting insertions to main chromosomes')
    main_chromosomes = map(str, range(1, 19+1)) + ['X', 'Y']
    insertions = insertions.ix[insertions['chromosome'].isin(main_chromosomes)]

    logging.info('Removing insertions with ulp < %d' % min_ulp)
    insertions = insertions.ix[insertions['unique_lp'] >= min_ulp]

    # Assign insertions ids and map to samples
    logging.info('Mapping insertions to samples')
    insertions['id'] = ['INS_%d' % ind for ind in  range(1, len(insertions)+1)]
    insertions = map_insertions_to_samples(insertions, barcode_sample_file, drop_barcode=True)

    logging.info('Done!')

    return insertions


def filter_contaminated_reads(reads, contaminant_file):
    filtered_reads = reads

    for contaminant in read_fasta(contaminant_file):
        filtered_reads = [r for r in filtered_reads if contaminant.seq not in r.seq]

    return filtered_reads
