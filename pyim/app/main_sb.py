

import HTSeq
import argparse
from os import path

from pyim.mapping.blat import insertion_statistics, blat_sequences, map_to_insertions
from pyim.protocols.sb.alignment import sb_alignments, genomic_sequences
from pyim.protocols.sb.stats import hist_alignment_quality

from pyim.plot.ggplot import ggplot_save
from pyim.plot.hist import hist_alignment_type, hist_barcode_alignment


PROJECT_DIR = path.join(path.dirname(path.abspath(__file__)), '../../')
DATA_DIR = path.join(PROJECT_DIR, 'data')


def _parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--fasta_file', required=True)
    parser.add_argument('--vector_file', default=path.join(DATA_DIR, 'vec.fa'))
    parser.add_argument('--barcode_file', default=path.join(DATA_DIR, 'SBbarcodes.fa'))
    parser.add_argument('--reference', required=True)
    parser.add_argument('--workdir', default='.')
    return parser.parse_args()


def main(fasta_file, vector_file, barcode_file, reference, workdir):
    # Load reads
    fasta_seqs = list(HTSeq.FastaReader(fasta_file))

    # Do alignments
    full_alignment, extra_info = sb_alignments(fasta_seqs, vector_file, barcode_file)

    #plot, _ = hist_alignment_quality(full_alignment, extra_info)
    #ggplot_save(plot, filepath=path.join(workdir, 'alignment_status.pdf'))

    # Check barcode alignments
    #plot, _ = _alignment_type_stats(full_alignment, extra_info, fasta_seqs, workdir)

    # Extract genomic sequences and do BLAT
    _, (bc_alignments, _) = extra_info
    genomic_seqs = genomic_sequences(full_alignment, bc_alignments)
    aligned_psl = blat_sequences(genomic_seqs, reference)

    insertions = map_to_insertions(aligned_psl, bc_alignments)



def _alignment_type_stats(full_alignment, extra_info, fasta_seqs, workdir):
    (vec_alignments, vec_unmapped), (bc_alignments, bc_unmapped) = extra_info

    # Check vector alignments
    hist_path = path.join(workdir, 'hist_aln_{vec}.pdf')

    plot, _ = hist_alignment_type(vec_alignments['T7'], unmapped=vec_unmapped['T7'], title='T7 alignment types')
    ggplot_save(plot, filepath=hist_path.format(vec='T7'))

    plot, _ = hist_alignment_type(vec_alignments['SB'], unmapped=vec_unmapped['SB'], title='SB alignment types')
    ggplot_save(plot, filepath=hist_path.format(vec='SB'))

    # Check barcode alignments
    plot, _ = hist_barcode_alignment(bc_alignments, bc_unmapped=bc_unmapped, title='Barcode alignments')
    ggplot_save(plot, filepath=path.join(workdir, 'hist_bc.pdf'))





if __name__ == '__main__':
    args = _parse_args()
    main(args.fasta_file, args.vector_file, args.barcode_file, args.workdir)
