import logging

import numpy

from pyim.io import read_fasta_filtered
from pyim.pipeline.shared.alignment import align_to_reference
from pyim.pipeline.sb.vectors import extract_genomic_seqs
from pyim.pipeline.sb.insertions import identify_insertions


logger = logging.getLogger('PyIM.SB')


def sb_pipeline(reads_file, run_name, reference, vector_aligners,
                vector_file, barcode_file, barcode_sample_file, work_dir,
                contaminant_file=None, threads=1, min_sequence_len=15, min_ulp=2):

    logger.info('Reading input reads')
    reads = read_fasta_filtered(reads_file, contaminant_file)

    logger.info('Aligning vectors and extracting genomic sequences')
    genomic_seqs = extract_genomic_seqs(reads, vector_aligners, vector_file, barcode_file, work_dir)

    logger.info('Aligning genomic sequences to reference genome')
    alignments = align_to_reference(genomic_seqs, reference, min_sequence_len, threads, work_dir)

    logger.info('Identifying insertions from alignments')
    insertions = identify_insertions(alignments, barcode_sample_file, min_ulp)
    insertions['id'] = insertions['id'].map(lambda id_: '%s.%s' % (run_name, id_))

    assert(not numpy.any(insertions['lp'] < insertions['unique_lp']))

    return insertions
