from pyim.alignment.vector.aligners import SequenceAligner


class ChainedReadAligner(SequenceAligner):

    def __init__(self, aligners):
        super(ChainedReadAligner, self).__init__()
        self.aligners = aligners

    def _align(self, queries, vector):
        unmapped, alignments = queries, {}
        
        for aligner in self.aligners:
            alns = aligner.align(unmapped, vector)
            alignments.update(alns)

            unmapped = [q for q in unmapped if q.name not in alignments]

        return alignments
