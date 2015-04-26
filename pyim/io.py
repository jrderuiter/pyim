from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
from builtins import (ascii, bytes, chr, dict, filter, hex, input,
                      int, map, next, oct, open, pow, range, round,
                      str, super, zip)

import skbio


def read_fasta(file_path):
    with file_path.open('r') as file_:
        for seq in skbio.io.read(file_, format='fasta'):
            yield seq
