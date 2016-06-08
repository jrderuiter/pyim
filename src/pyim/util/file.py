import pysam
from functools import reduce


def _make_gen(reader):
    b = reader(1024 * 1024)
    while b:
        yield b
        b = reader(1024*1024)


def count_lines(file_path):
    f = open(file_path, 'rb')
    f_gen = _make_gen(f.raw.read)
    return sum(buf.count(b'\n') for buf in f_gen)


def count_fasta_entries(file_path):
    f = open(file_path, 'rb')
    f_gen = _make_gen(f.raw.read)
    return sum(buf.count(b'>') for buf in f_gen)


def count_bam_entries(file_path):
    # From Biostars at https://www.biostars.org/p/1890/.
    # Could be faster for sorted/index bam files using idxstats.
    reduce(lambda x, y: x + y,
           [eval('+'.join(l.rstrip('\n').split('\t')[2:]))
            for l in pysam.idxstats(file_path)])
