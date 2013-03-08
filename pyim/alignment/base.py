
import pandas as pd


class Alignment(object):

    def __init__(self, q_name, q_start, q_end, q_seq, t_name, t_start, t_end, t_seq, identity, score, alnStr=None):
        self.query = self._internal_dict(q_name, q_start, q_end, q_seq)
        self.target = self._internal_dict(t_name, t_start, t_end, t_seq)

        self.identity = identity
        self.score = score
        self.alignment = alnStr

    def _internal_dict(self, name, start, end, seq):
        return {
            'name': name,
            'start': start,
            'end': end,
            'seq': seq
        }

    def in_query(self):
        readSeq = list(self.query['seq'])
        start, end = self.query['start'], self.query['end']
        vecName = self.target['name']

        size = end - start
        paddingSize = size - len(vecName) - 3
        vecStr = '[-' + vecName + ('-' * paddingSize) + ']'

        readSeq[start:end] = list(vecStr)
        print(''.join(readSeq))

    def as_dict(self):
        dict_ = {}
        for key, value in self.query.items():
            dict_['query_' + key] = value
        for key, value in self.target.items():
            dict_['target_' + key] = value
        dict_['identity'] = self.identity
        dict_['score'] = self.score
        dict_['alignment'] = self.alignment
        return dict_

    @classmethod
    def as_dataframe(cls, alignments):
        if isinstance(alignments, dict):
            alignmentDicts = [al.as_dict() for al in alignments.itervalues()]
            alignmentFrame = pd.DataFrame(alignmentDicts, index=alignments.keys())
        else:
            alignmentDicts = [al.as_dict() for al in alignments]
            alignmentFrame = pd.DataFrame(alignmentDicts)
        alignmentFrame.reindex_axis([], axis=1)
        return alignmentFrame

