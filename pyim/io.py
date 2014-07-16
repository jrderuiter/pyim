from __future__ import print_function

import os
from os import path
from itertools import groupby
from collections import namedtuple
from pyim.util import chunks


Sequence = namedtuple('Sequence', ['name', 'seq'])


def read_fasta(file_path, as_dict=False):
    if as_dict:
        return { seq.name: seq for seq in fasta_generator(file_path) }
    else:
        return (seq for seq in fasta_generator(file_path))


def read_fasta_filtered(reads_file, contaminant_file=None):
    reads = list(read_fasta(reads_file))

    if contaminant_file is not None:
        for contaminant in read_fasta(contaminant_file):
            reads = [r for r in reads if contaminant.seq not in r.seq]

    return reads


def fasta_generator(file_path):
    with open(file_path, 'r') as file_:
        fa_iter = (x[1] for x in groupby(file_, lambda line: line[0] == ">"))
        for header in fa_iter:
            header = next(header)[1:].strip().split(' ')[0]     # drop the ">" and select first element
            seq = "".join(s.strip() for s in next(fa_iter))     # join all sequence lines to one.
            yield Sequence(header, seq)


def write_fasta(file_path, seqs, width=80):
    with open(file_path, 'w') as file_:
        for seq in seqs:
            print('>%s' % seq.name, file=file_)
            for chunk in chunks(seq.seq, width):
                print(chunk, file=file_)
    return file_path


def write_insertions_to_gff(insertions, file_path):
    with open(file_path, 'w') as gff_file:
        gff_str = '{chrom}\t.\tinsertion\t{start}\t{end}\t.\t{strand}\t.\t{info}'
        info_str = 'ID={name}({id});NAME={id};LP={lp};uLP={ulp}'

        for i, (_, row) in enumerate(insertions.iterrows()):
            name = row['ins_id'] if 'ins_id' in row else 'INS_' + str(i)

            info_str_fmt = info_str.format(name=name, id=row['sample'], lp=row['lp'], ulp=row['unique_lp'])
            gff_str_fmt = gff_str.format(chrom=row['chromosome'],  strand=row['strand'], info=info_str_fmt,
                                     start=row['location'] - 10, end=row['location'] + 10)

            print(gff_str_fmt, file=gff_file)

    return file_path


def makedirs_safe(dir_path):
    if not path.exists(dir_path):
        os.makedirs(dir_path)
    return dir_path


def write_frame(frame, file_path, sep='\t', index=False):
    file_dir = os.path.dirname(file_path)
    if not os.path.exists(file_dir):
        os.makedirs(file_dir)
    frame.to_csv(file_path, sep=sep, index=index)