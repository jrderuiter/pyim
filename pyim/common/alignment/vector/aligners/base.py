from pyim.common.alignment.vector.model import VectorAlignment


class VectorAligner(object):
    """docstring for VectorAligner"""

    def __init__(self, seq_aligner=None):
        super(VectorAligner, self).__init__()
    
        if seq_aligner is None:
            seq_aligner = ExactSequenceAligner()

        self.aligner = seq_aligner
        
    def align(self, sequences, vector):
        return self.aligner.align(sequences, vector)

    def align_multiple(self, sequences, vectors):
        return self.aligner.align_multiple(sequences, vectors)
                                 

class SequenceAligner(object):
    """docstring for SequenceAligner"""

    def __init__(self):
        super(SequenceAligner, self).__init__()
        
    def align(self, queries, vector):
        return self._align(queries, vector)

    def _align(self, queries, vector):
        raise NotImplementedError

    def align_multiple(self, queries, vectors):
        alignments = self._align_multiple(queries, vectors)
        merged = self._merge_alignments(alignments)
        return merged

    def _align_multiple(self, queries, vectors):
        return [self.align(queries, v) for v in vectors]

    def _merge_alignments(self, alignments):
        merged, duplicates = {}, {}
        
        # Merge all alignments, marking duplicates. Note we
        # don't collect the first entry of a duplicate (which
        # is already inserted into merged) until the next step,
        # as this allows us to avoid an extra key lookup in
        # the duplicates dictionary.
        for vec_alns in alignments:
            for aln in vec_alns.values():
                query_name = aln.q_name

                if query_name not in merged:
                    merged[query_name] = aln
                else:
                    duplicates = [aln]
        
        # Now we collect the first entry of a duplicate.
        for query_name in duplicates.keys():
            duplicates[query_name].append(merged[query_name])

        # And finally we attempt to resolve the duplicates.
        # If we can't resolve the duplicates (resolve function
        # returns None), then we delete the entire alignment.
        for query_name, dupl in duplicates.items():
            resolved = self._resolve_duplicates(dupl)

            if resolved is None:
                print("WARNING: Can't resolve multiple alignments " + \
                      "for %s, dropping alignments." % query_name)
                del merged[query_name]
            else:
                merged[query_name] = resolved
        
        return merged

    def _resolve_duplicates(self, duplicates):
        return None


class ExactSequenceAligner(SequenceAligner):
    """docstring for ExactSequenceAligner"""

    def _align(self, queries, vector):
        v_seq = vector.sequence
        v_len = len(v_seq)
        
        alignments = {}
        for q in queries:
            try:
                ind = q.sequence.index(v_seq)
                aln = VectorAlignment(
                    q_name=q.name, q_start=ind, q_end=ind + v_len,
                    v_name=vector.name, v_start=0, v_end=v_len)
                alignments[q.name] = aln
            except ValueError:
                pass
            
        return alignments
