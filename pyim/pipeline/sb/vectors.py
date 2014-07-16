from os import path

import brewer2mpl
import pandas
import numpy
from matplotlib import pyplot as plt

from pyim.io import read_fasta, makedirs_safe, write_frame

## Set style to make plots slightly more pretty.
pandas.options.display.mpl_style = 'default'


def extract_genomic_seqs(reads, aligners, vector_file, barcode_file, work_dir):
    """Extract genomic sequences from SB reads.

    Uses the supplied aligners to recognise the SB, T7 and barcode sequences
    in the supplied set of sequence reads. If all three sequences are identified
    the genomic sequence between the SB and the T7 sequence is extracted from
    the original read sequence. This genomic sequence is returned in a frame
    together with the original read id and the identified sample barcode.

    Statistics and intermediate results are written to 'work_dir/_vectors'.
    These files may be used for debugging purposes.

    :Parameters:
        - `reads` (Sequence list)   - SB input reads from the 454 sequencer.
        - `aligners` (Aligner dict) - Dictionary containing the vector aligners that should be used.
                                      The dictionary should contain aligners for 't7', 'sb' and 'bc'.
        - `vector_file` (string)    - Path to the file containing the SB and T7 vectors.
        - `barcode_file` (string)   - Path to the file containing the barcode vectors.
        - `work_dir` (string)       - Working directory to write files to.

    :Returns:
        Frame containing the columns query_name, bc_name and genomic_seq. Each row contains
        the genomic sequence for the read identified by query_name. Each entry is tagged with
        the barcode that was identified from the read.

    :Returns Type:
        pandas.DataFrame

    :Examples:

    """

    # Create the working directory
    work_dir = makedirs_safe(path.join(work_dir, '_vectors'))

    # Perform vector alignments and write/plot statistics
    vector_alignments = align_vectors(reads, aligners, vector_file, barcode_file)
    write_frame(vector_alignments, path.join(work_dir, 'vector_alignments.txt'))
    _plot_alignment_stats(vector_alignments, work_dir)

    # Extract genomic sequences
    genomic_seqs = alignments_to_genomic_seqs(vector_alignments)
    write_frame(genomic_seqs, path.join(work_dir, 'genomic_sequences.txt'))
    _plot_sequence_hist(genomic_seqs, work_dir)

    return genomic_seqs


def align_vectors(reads, aligners, vector_file, barcode_file):
    # Read vectors to which we will align our reads
    vector_seqs = read_fasta(vector_file, as_dict=True)
    barcode_seqs = list(read_fasta(barcode_file, as_dict=False))

    # Perform the actual alignment to each of the sequences
    sb_alignment, sb_unaligned = align_vector(reads, vector_seqs['SB'], aligners['sb'])
    t7_alignment, t7_unaligned = align_vector(reads, vector_seqs['T7'], aligners['t7'])
    bc_alignment, bc_unaligned = align_barcodes(reads, barcode_seqs, aligners['bc'])

    # Extract the genomic sequences and bar-codes using these alignments
    return merge_alignments(sb_alignment, t7_alignment, bc_alignment)


def align_vector(reads, seq, aligner):
    return aligner.align_target(reads, seq)  # Returns alignment, unmapped


def align_barcodes(reads, barcode_seqs, aligner=None):
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
    merged = sb_sub.merge(t7_sub, on='query_name', how='outer')
    merged = merged.merge(bc_sub, on='query_name', how='outer')

    return merged


def alignments_to_genomic_seqs(merged_alignment):

    # Only use complete entries for extracting genomic sequences
    merged_alignment = merged_alignment.dropna(axis=0, how='any')

    # Generate list of genomic sequences
    items = zip(merged_alignment['query_seq'],
                merged_alignment['sb_end'].astype(numpy.int64),
                merged_alignment['t7_start'].astype(numpy.int64))
    seqs = [seq[sb_end:t7_start] for seq, sb_end, t7_start in items]

    # Return genomic sequences as data frame
    frame = {'query_name': merged_alignment['query_name'],
             'bc_name': merged_alignment['bc_name'],
             'genomic_seq': seqs}
    frame = pandas.DataFrame.from_dict(frame)

    return frame[['query_name', 'bc_name', 'genomic_seq']]


def _plot_alignment_stats(merged_alignment, work_dir):
    work_dir = makedirs_safe(path.join(work_dir, 'plots'))

    type_cols = ['sb_alignment_type', 't7_alignment_type', 'bc_alignment_type']
    aln_types = merged_alignment[type_cols].fillna('none')

    base_path = path.join(work_dir, '%s_aln_hist.pdf')
    _plot_aln_hist(aln_types, 't7_alignment_type', base_path % 't7', title='T7 alignments')
    _plot_aln_hist(aln_types, 'sb_alignment_type', base_path % 'sb', title='SB alignments')
    _plot_aln_hist(aln_types, 'bc_alignment_type', base_path % 'bc', title='Barcode alignments')

    _plot_bc_hist(merged_alignment, path.join(work_dir, 'bc_counts.pdf'))


def _plot_aln_hist(alignments, column, file_path, title=''):
    counts = alignments[column].value_counts()
    colors = brewer2mpl.get_map('Set2', 'qualitative', 8).mpl_colors
    _bar_plot(counts, file_path, title, fig_size=(6, 6), color=colors, margin_bottom=0.35)


def _plot_bc_hist(alignments, file_path, title='Barcode counts'):
    counts = alignments['bc_name'].value_counts(sort=False)
    counts = counts[counts.index.order()]
    _bar_plot(counts, file_path, title, fig_size=(10, 6))


def _bar_plot(counts, file_path, title='', fig_size=None, color=None, margin_bottom=None):
    plt.ioff()
    fig = plt.figure(figsize=fig_size)

    ax = counts.plot(kind='bar', title=title, color=color)
    ax.set_xlim((-0.5, len(counts)-0.5))

    if margin_bottom is not None:
        fig.subplots_adjust(bottom=margin_bottom)

    fig.savefig(file_path)
    plt.close(fig)


def _plot_sequence_hist(genomic_seqs, work_dir):
    ## Prepare path and inputs
    work_dir = makedirs_safe(path.join(work_dir, 'plots'))
    file_path = path.join(work_dir, 'genomic_seq_length.pdf')
    seq_lens = genomic_seqs['genomic_seq'].str.len()

    ## Plot figure
    plt.ioff()
    fig = plt.figure(figsize=(10, 6))

    bin_size = 10
    ax = seq_lens.hist(bins=range(0, seq_lens.max() + bin_size, bin_size))
    ax.set_title('Genomic sequence length distribution')

    fig.savefig(file_path)
    plt.close(fig)
