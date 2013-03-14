
import pandas as pd
from pyim.alignment.model import Alignment


class ReadAligner(object):

    def align_target(self, fasta_seqs, target_seq):
        raise NotImplementedError

    def align_targets(self, fasta_seqs, target_seqs):
        target_alignments = [self.align_target(fasta_seqs, tgt) for tgt in target_seqs]
        return pd.concat(target_alignments)

    def _to_frame(self, alignments):
        if len(alignments) == 0:
            return None

        if isinstance(alignments[0], Alignment):
            alignments = [aln.__dict__ for aln in alignments]

        frame = pd.DataFrame(alignments)
        frame.set_index('query_name', inplace=True)

        return frame
