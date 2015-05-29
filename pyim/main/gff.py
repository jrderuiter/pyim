from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
from builtins import (ascii, bytes, chr, dict, filter, hex, input,
                      int, map, next, oct, open, pow, range, round,
                      str, super, zip)
from future.utils import native_str

from argparse import ArgumentParser
from pathlib import Path

import pandas as pd


def setup_parser():
    parser = ArgumentParser(prog='pyim-gff')

    parser.add_argument('insertions', type=Path)
    parser.add_argument('output_gff', type=Path)

    return parser


def _ins_to_gff(ins, size=1000):
    assert isinstance(ins.strand, int)

    attrs = [i for i in ins.index if i not in
             {'insertion_id', 'seqname', 'location', 'strand'}]

    attr_dict = {attr: ins[attr] for attr in attrs}

    attr_dict['id'] = ins.insertion_id
    attr_dict['name'] = ins.insertion_id

    attr_keys = sorted(attr_dict.keys())
    attr_str = ';'.join(('{} {}'.format(k, attr_dict[k]) for k in attr_keys))

    return {
        'seqname': ins.seqname,
        'source': '.',
        'feature': 'insertion',
        'start': int(ins.location - (size / 2)),
        'end': int(ins.location + (size / 2)),
        'score': '.',
        'strand': '+' if ins.strand == 1 else '-',
        'frame': '.',
        'attribute': attr_str
    }


def main():
    parser = setup_parser()
    args = parser.parse_args()

    # Read input.
    ins_frame = pd.read_csv(str(args.insertions), sep=native_str('\t'),
                            dtype={'seqname': str, 'location': int})

    # Transform to gff frame.
    gff_frame = pd.DataFrame.from_records(
        (_ins_to_gff(r) for _, r in ins_frame.iterrows()),
        columns=['seqname', 'source', 'feature', 'start', 'end',
                 'score', 'strand', 'frame', 'attribute'])
    gff_frame = gff_frame.sort(['seqname', 'start', 'end'])

    # Write output.
    gff_frame.to_csv(str(args.output_gff), sep=native_str('\t'),
                     index=False, header=False)


if __name__ == '__main__':
    main()
