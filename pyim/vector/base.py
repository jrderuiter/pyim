
import pandas
from collections import namedtuple
from pyim.io import read_fasta
from pyim.alignment.aligners.base import ExactReadAligner


class VectorAligner(object):

    def align(self, reads, vector_file, barcode_file):
        pass

    def _read_seqs(self, fasta_file, name_mask=None):
        seqs = { seq.name: seq for seq in read_fasta(fasta_file) }
        if name_mask is not None:
            seqs =  { k: v for k, v in self._read_seqs(fasta_file).items() if k in name_mask }
        return seqs

    def _align_vector(self, reads, vector, aligner=None):
        if aligner is None: aligner = ExactReadAligner()
        alignment, unmapped = aligner.align_target(reads, vector)
        return alignment, unmapped

    def _align_barcodes(self, reads, barcodes, aligner=None):
        if aligner is None: aligner = ExactReadAligner()

        alignments, _ = aligner.align_targets(reads, barcodes)
        alignments = pandas.concat(alignments.values(), ignore_index=True)

        is_mapped = pandas.Series([read.name for read in reads]).isin(alignments['query_name'])
        unmapped = [reads[i] for i, mapped in enumerate(is_mapped) if not mapped]

        return alignments, unmapped


AlignmentResult = namedtuple('AlignmentResult', ['genomic_sequences', 'combined_alignments',
                                                 'vector_alignments', 'vector_unmapped',
                                                 'barcode_alignments', 'barcode_unmapped'])

#class AlignmentResult(_result_tup):
#    def __new__(cls, ):
#        if extra is None: extra = {}
#        return super(AlignmentResult, cls).__new__(cls, genomic_seqs, vec_alignments, bc_alignments, extra)

