import os
from itertools import groupby
from os import path

from pyim.common.util import chunks
from pyim.common.model import Sequence


def makedirs_safe(dir_path):
    if not path.exists(dir_path):
        os.makedirs(dir_path)
    return dir_path


def read_fasta(file_path, as_dict=False):
    if as_dict:
        return { seq.name: seq for seq in fasta_generator(file_path) }
    else:
        return (seq for seq in fasta_generator(file_path))


def fasta_generator(file_path):
    with open(file_path, 'r') as file_:
        fa_iter = (x[1] for x in groupby(file_, lambda line: line[0] == ">"))
        for header in fa_iter:
            # Drop the ">" and select first element.
            header = next(header)[1:].strip()
            header = header.split(' ')[0]

            # Join all sequence lines to one.
            seq = "".join(s.strip() for s in next(fa_iter))
            
            yield Sequence(name=header, sequence=seq)


def write_fasta(sequences, file_path, width=80):
    with open(file_path, 'w') as file_:
        for seq in sequences:
            print('>%s' % seq.name, file=file_)
            for chunk in chunks(seq.sequence, width):
                print(chunk, file=file_)
    return file_path
