from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
from builtins import (ascii, bytes, chr, dict, filter, hex, input,
                      int, map, next, oct, open, pow, range, round,
                      str, super, zip)

from itertools import chain, repeat

import numpy as np

from skbio import BiologicalSequence
from skbio.io import register_reader, register_writer


@register_reader('faster_fastq')
def _fastq_to_generator(fh, constructor=BiologicalSequence,
                        qual_in_description=True, phred_offset=33):
    for id_, seq, _, qual in grouper(4, fh):
        if qual_in_description:
            yield constructor(id=id_.strip(), sequence=seq.strip(),
                              description=qual.strip())
        else:
            yield constructor(id=id_.strip(), sequence=seq.strip(),
                              quality=decode_qual(qual.strip(), phred_offset))


def decode_qual(qual_str, phred_offset=33):
    return np.fromstring(qual_str, dtype='uint8') - phred_offset


@register_writer('faster_fastq')
def _generator_to_fastq(obj, fh, qual_in_description=True):

    for sequence in obj:
            quality = sequence.description if qual_in_description \
                else sequence.quality

            print(sequence.id, file=fh)
            print(sequence.sequence, file=fh)
            print('+', file=fh)
            print(quality, file=fh)


@register_writer('faster_fastq', BiologicalSequence)
def _biological_sequence_to_fastq(obj, fh, **kwargs):
    _sequences_to_fastq([obj], fh, **kwargs)


def _sequences_to_fastq(obj, fh, **kwargs):
    def seq_gen():
        for seq in obj:
            yield seq

    _generator_to_fastq(seq_gen(), fh, **kwargs)


def grouper(n, iterable, pad_value=None):
    """ Group iterable into chunks of n.

    grouper(3, 'abcdefg', 'x') --> ('a','b','c'), ('d','e','f'), ('g','x','x')"

    :param n:
    :param iterable:
    :param pad_value:
    :return:
    """
    return zip(*[chain(iterable, repeat(pad_value, n-1))]*n)
