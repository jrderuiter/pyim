
import HTSeq
import pandas

from pyim.alignment.aligners.base import ReadAligner


class ChainedReadAligner(ReadAligner):

    def __init__(self, aligners):
        super(ChainedReadAligner, self).__init__()
        self.aligners = aligners

    def align_target(self, fasta_seqs, target_seq):
        unmapped_reads = fasta_seqs

        alignments = []
        for aligner in self.aligners:
            aln_alignments, unmapped_reads = aligner.align_target(unmapped_reads, target_seq)
            alignments.append(aln_alignments)

        return pandas.concat(alignments, ignore_index=True), unmapped_reads


class TruncatedTargetAligner(ReadAligner):

    def __init__(self, aligner, trunc_length):
        super(TruncatedTargetAligner, self).__init__()
        self.aligner = aligner
        self.trunc_length = trunc_length

    def align_target(self, fasta_seqs, target_seq):
        target_len = len(target_seq)

        trunc_target = HTSeq.Sequence(target_seq.seq[:self.trunc_length*-1], target_seq.name)
        alignments, unmapped_reads = self.aligner.align_target(fasta_seqs, trunc_target)

        # Fix some fields to indicate truncation has taken place
        # TODO: adjust cigar?
        type_suffix = "_truncated[%d]" % self.trunc_length
        alignments['type'] = alignments.iloc[0]['type'] + type_suffix
        alignments['identity'] = (alignments['identity'] * (target_len - self.trunc_length)) / target_len
        alignments['score'] = alignments['identity'].mul(100).round()

        return alignments, unmapped_reads