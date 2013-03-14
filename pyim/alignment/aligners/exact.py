
from pyim.alignment.aligners.base import ReadAligner
from pyim.alignment.model import Alignment


class ExactReadAligner(ReadAligner):

    def align_target(self, fasta_seqs, target):
        tgt_seq, tgt_len = target.seq, len(target)
        hits = [s for s in fasta_seqs if tgt_seq in s.seq]

        alignments = []

        for hit in hits:
            start = hit.seq.index(tgt_seq)
            cigar_str = '%dM' % tgt_len
            alignment = Alignment(hit.seqId, start, start + tgt_len, hit.seq,
                                  target.seqId, 0, tgt_len, tgt_seq,
                                  100, 1.0, cigar_str, 'exact')
            alignments.append(alignment)

        return self._to_frame(alignments)

