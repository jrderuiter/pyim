
import itertools


def read_fasta(filePath):
    with open(filePath, 'r') as fastaFile:
        # Ditch boolean as we know groups alternate
        sectIter = (x[1] for x in itertools.groupby(fastaFile, _is_header))
        for header in sectIter:
            seqId, metadata = _parse_header(next(header).strip())
            seq = ''.join([s.strip() for s in next(sectIter)])
            yield FastaSequence(seqId, seq, metadata)


def _parse_header(header):
    splitHeader = header.split(' ')
    seqId = splitHeader[0][1:]

    metadata = {}
    for entry in splitHeader[1:]:
        name, value = entry.split('=')
        metadata[name] = value

    return seqId, metadata


def _is_header(line):
    return line[0] == '>'


class FastaSequence(object):

    def __init__(self, seqId, seq, metadata=None):
        self.seqId = seqId
        self.seq = seq
        self.metadata = metadata

    def __len__(self):
        return len(self.seq)

    def __str__(self):
        return '<Fasta Sequence %s: %s>' % (self.seqId, self.seq)

    def __repr__(self):
        return self.__str__()

