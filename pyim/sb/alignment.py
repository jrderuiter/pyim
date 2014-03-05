import pandas

from pyim.alignment.vector.aligners import ExactReadAligner
from pyim.io import Sequence, read_fasta

DEFAULT_ALIGNER = ExactReadAligner


def align(reads, vector_file, barcode_file, sb_aligner=None, t7_aligner=None, bc_aligner=None):
    # Read vectors to which we will align our reads
    vector_seqs = read_fasta(vector_file, as_dict=True)
    barcode_seqs = list(read_fasta(barcode_file, as_dict=False))

    # Perform the actual alignment to each of the sequences
    sb_alignment, sb_unaligned = align_sb(reads, vector_seqs['SB'], sb_aligner)
    t7_alignment, t7_unaligned = align_t7(reads, vector_seqs['T7'], t7_aligner)
    bc_alignment, bc_unaligned = align_barcodes(reads, barcode_seqs, bc_aligner)

    # Extract the genomic sequences and bar-codes using these alignments
    merged_alignment = merge_alignments(sb_alignment, t7_alignment, bc_alignment)
    genomic_seqs, barcode_mapping = extract_genomic_sequences(merged_alignment)

    return genomic_seqs, barcode_mapping


def align_sb(reads, sb_seq, aligner):
    return _align_vector(reads, sb_seq, aligner)


def align_t7(reads, t7_seq, aligner):
    return _align_vector(reads, t7_seq, aligner)


def _align_vector(reads, seq, aligner):
    if aligner is None:
        aligner = DEFAULT_ALIGNER()
    return aligner.align_target(reads, seq)  # Returns alignment, unmapped


def align_barcodes(reads, barcode_seqs, aligner=None):
    if aligner is None:
        aligner = DEFAULT_ALIGNER()

    alignments, _ = aligner.align_targets(reads, barcode_seqs)
    alignments = pandas.concat(alignments.values(), ignore_index=True)

    is_mapped = pandas.Series([read.name for read in reads]).isin(alignments['query_name'])
    unmapped = [reads[i] for i, mapped in enumerate(is_mapped) if not mapped]

    return alignments, unmapped


def merge_alignments(sb_alignment, t7_alignment, bc_alignment):
    sel_columns = ['query_name', 'query_start', 'query_end', 'type']

    # Select columns and rename before merging
    sb_sub = sb_alignment[sel_columns]
    sb_sub.columns = ['query_name', 'sb_start', 'sb_end', 'sb_alignment_type']

    t7_sub = t7_alignment[sel_columns]
    t7_sub.columns = ['query_name', 't7_start', 't7_end', 't7_alignment_type']

    bc_sub = bc_alignment[sel_columns + ['target_name', 'query_seq']]
    bc_sub.columns = ['query_name', 'bc_start', 'bc_end', 'bc_alignment_type', 'bc_name', 'query_seq']

    # Merge frames into main alignment frame, use inner join to only
    # select those reads that have all three alignments.
    merged = sb_sub.merge(t7_sub, on='query_name', how='inner')
    merged = merged.merge(bc_sub, on='query_name', how='inner')

    return merged


def extract_genomic_sequences(merged_alignment):
    # Fetch Series objects for iteration
    query_names = merged_alignment['query_name']
    query_seqs = merged_alignment['query_seq']
    sb_ends = merged_alignment['sb_end']
    t7_starts = merged_alignment['t7_start']

    # Collect sequences
    seqs = []
    for i, query_seq in enumerate(query_seqs):
        seq_str = query_seq[sb_ends.iloc[i]:t7_starts.iloc[i]]
        seq = Sequence(query_names.iloc[i], seq_str)
        seqs.append(seq)

    # Return sequences together with mapping of sequence to barcode
    return seqs, dict(zip(query_names, merged_alignment['bc_name']))


