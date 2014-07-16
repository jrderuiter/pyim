from os import path

from pyim.alignment.genomic.bowtie2 import Bowtie2Aligner
from pyim.io import makedirs_safe, Sequence, write_frame


def align_to_reference(genomic_seqs, reference, min_sequence_len, threads, work_dir):
    """Aligns genomic sequences to the supplied reference genome.

    :Parameters:
        - `genomic_seqs` (pandas.DataFrame) - DataFrame containing the genomic sequences. Frame contains
                                              at least the two columns 'query_name' and 'genomic_seq'.
        - `reference` (string)              - Path to the fasta genome reference sequence.
        - `min_sequence_len` (int)          - Minimum length for a sequence to be included for alignment.
        - `threads` (int)                   - Number of threads to be used in aligning sequences.
        - `work_dir` (string)               - Working directory to write files to.

    :Returns:
        DataFrame containing the sequence alignments, merged with the input sequence frame.
        Note that genomic sequences without an alignment are dropped from the output.

    :Returns Type:
        pandas.DataFrame

    :Examples:

    """

    work_dir = makedirs_safe(path.join(work_dir, '_alignment'))

    genomic_seqs = genomic_seqs.ix[genomic_seqs['genomic_seq'].str.len() >= min_sequence_len]
    seq_objs = _frame_to_seqs(genomic_seqs)

    genomic_aligner = Bowtie2Aligner(work_dir, num_cores=threads)
    hits, _, _, _ = genomic_aligner.align(seq_objs, reference)

    merged = hits.merge(genomic_seqs, on='query_name', how='inner')
    write_frame(merged, path.join(work_dir, 'alignments.txt'))

    return merged


def _frame_to_seqs(genomic_seqs):
    """Converts DataFrame of genomic sequences to list of Sequence objects.

    :Parameters:
        - `genomic_seqs` (pandas.DataFrame) - DataFrame containing the genomic sequences. Frame contains
                                              at least the two columns 'query_name' and 'genomic_seq'.

    :Returns:
        List of Sequence objects reflecting the sequences of the input.

    :Returns Type:
        List of Sequences.

    :Examples:

    """

    items = zip(genomic_seqs['query_name'], genomic_seqs['genomic_seq'])
    return [Sequence(name, seq) for name, seq in items]