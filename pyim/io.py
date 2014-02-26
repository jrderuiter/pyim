
from __future__ import print_function

from itertools import groupby
from collections import namedtuple
from pyim.util import chunks


Sequence = namedtuple('Sequence', ['name', 'seq'])

def read_fasta(file_path):
    with open(file_path, 'r') as file_:
        fa_iter = (x[1] for x in groupby(file_, lambda line: line[0] == ">"))
        for header in fa_iter:
            header = header.next()[1:].strip().split(' ')[0]             # drop the ">" and select first element
            seq = "".join(s.strip() for s in fa_iter.next())     # join all sequence lines to one.
            yield Sequence(header, seq)


def write_fasta(file_path, seqs, width=80):
    with open(file_path, 'w') as file_:
        for seq in seqs:
            print('>%s' % seq.name, file=file_)
            for chunk in chunks(seq.seq, width):
                print(chunk, file=file_)
    return file_path

