from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
from builtins import (ascii, bytes, chr, dict, filter, hex, input,
                      int, map, next, oct, open, pow, range, round,
                      str, super, zip)
from future.utils import native_str

from argparse import ArgumentParser
from pathlib import Path

import pandas as pd

from pyim.cis.cimpl import cimpl, get_cis, get_cis_mapping


def setup_parser():
    parser = ArgumentParser(prog='pyim-cis')

    parser.add_argument('input', type=Path)
    parser.add_argument('output', type=Path)

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--pattern', default=None)
    group.add_argument('--system', choices={'SB'}, default=None)

    parser.add_argument('--genome', choices={'mm10'}, default='mm10')
    parser.add_argument('--chromosomes', nargs='+', default=None)
    parser.add_argument('--scales', nargs='+', type=int, default=30000)

    parser.add_argument('--alpha', type=float, default=0.05)
    parser.add_argument('--iterations', type=int, default=1000)
    parser.add_argument('--threads', type=int, default=1)
    parser.add_argument('--verbose', default=False, action='store_true')

    return parser


def main():
    parser = setup_parser()
    args = parser.parse_args()

    # Read frame.
    ins_frame = pd.read_csv(str(args.input), sep=native_str('\t'))

    # Run cimpl.
    cimpl_obj = cimpl(ins_frame, scales=args.scales, genome=args.genome,
                      system=args.system, pattern=args.pattern,
                      chromosomes=args.chromosomes, iterations=args.iterations,
                      threads=args.threads, verbose=args.verbose)

    # Extract cis and cis mapping from object.
    cis = get_cis(cimpl_obj, alpha=args.alpha, mul_test=True)
    cis_mapping = get_cis_mapping(cimpl_obj, cis_frame=cis)

    ins_annotated = pd.merge(ins_frame, cis_mapping, on='insertion_id')

    # Write out outputs.
    cis.to_csv(str(args.output.with_suffix('.cis.txt')),
               sep=native_str('\t'), index=False)
    ins_annotated.to_csv(str(args.output), sep=native_str('\t'), index=False)


if __name__ == '__main__':
    main()
