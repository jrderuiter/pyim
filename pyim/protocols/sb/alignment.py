

import HTSeq
import pandas

from pyim.alignment.aligners.base import ExactReadAligner
from pyim.alignment.aligners.composed import TruncatedTargetAligner, ChainedReadAligner


def sb_alignments(fasta_seqs, vector_file, barcode_file):
    sb_vectors = [seq for seq in HTSeq.FastaReader(vector_file) if seq.name in ['T7', 'SB']]
    vec_alignments, vec_unmapped = sb_vector_alignments(fasta_seqs, sb_vectors)

    bc_vectors = list(HTSeq.FastaReader(barcode_file))
    bc_alignments, bc_unmapped = sb_barcode_alignments(fasta_seqs, bc_vectors)

    merged = combine_alignments(vec_alignments, bc_alignments)
    extra = ((vec_alignments, vec_unmapped), (bc_alignments, bc_unmapped))

    return merged, extra


def sb_vector_alignments(fasta_seqs, sb_vectors):
    trunc_aligner = TruncatedTargetAligner(ExactReadAligner(), 1)
    comp_aligner = ChainedReadAligner([ExactReadAligner(), trunc_aligner])

    alignments, unmapped_seqs = comp_aligner.align_targets(fasta_seqs, sb_vectors)

    return alignments, unmapped_seqs


def sb_barcode_alignments(fasta_seqs, bc_vectors):
    aligner =  ExactReadAligner()
    alignments, unmapped_seqs = aligner.align_targets(fasta_seqs, bc_vectors)

    # Return only reads that never aligned anywhere!
    bc_aligned_queries = pandas.Index(pandas.concat(alignments, ignore_index=True)['query_name'])
    bc_unmapped_seqs = [s for s in fasta_seqs if s.name not in bc_aligned_queries]

    return alignments, bc_unmapped_seqs


def combine_alignments(vec_alignments, bc_alignments):
    bc_alignments = pandas.concat(bc_alignments.values(), ignore_index=True)

    sel_columns = ['query_name', 'query_start', 'query_end', 'type']
    merged = vec_alignments['SB'][sel_columns].merge(vec_alignments['T7'][sel_columns],
                                                     on='query_name', how='inner', suffixes=['_SB', '_T7'])
    merged = merged.merge(bc_alignments[sel_columns + ['target_name', 'query_seq']], on='query_name', how='inner')

    return merged


def genomic_sequences(vector_barcode_frame, minlength=0):
    query_names = vector_barcode_frame['query_name']
    query_seqs = vector_barcode_frame['query_seq']
    sb_ends = vector_barcode_frame['query_end_SB']
    t7_starts = vector_barcode_frame['query_start_T7']

    genomic_seqs = []
    for i, query_seq in enumerate(query_seqs):
        genomic_seq = query_seq[sb_ends.iloc[i]:t7_starts.iloc[i]]
        if len(genomic_seq) >= minlength:
            genomic_seqs.append(HTSeq.Sequence(genomic_seq, query_names.iloc[i]))

    return genomic_seqs
