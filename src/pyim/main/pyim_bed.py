import argparse
from collections import OrderedDict
from pathlib import Path

import numpy as np
import pandas as pd

from pyim.model import Insertion


def main():
    """Main function for pyim-cis."""

    args = parse_args()

    # Read insertions.
    ins_frame = Insertion.from_csv(args.insertions, sep='\t', as_frame=True)

    # Drop any columns if needed.
    if args.drop is not None:
        ins_frame = ins_frame.drop(args.drop, axis=1)
        ins_frame = ins_frame.drop_duplicates()

    # Convert to BED frame.
    start = (ins_frame['position'] - (args.width // 2)).astype(int)
    end = (ins_frame['position'] + (args.width // 2)).astype(int)
    strand = ins_frame['strand'].map({1: '+', -1: '-', np.nan: '.'})
    color = strand.map({'+': '0,0,255', '-': '255,0,0', '.': '60,60,60'})

    bed_frame = pd.DataFrame(
        OrderedDict([
            ('chrom', ins_frame['chromosome']),
            ('chromStart', start),
            ('chromEnd', end),
            ('name', ins_frame['id']),
            ('score', ins_frame['support']),
            ('strand', strand),
            ('thickStart', start),
            ('thickEnd', end),
            ('itemRgb', color)
        ])
    )  # yapf: disable

    # Write output.
    bed_frame.to_csv(str(args.output), sep='\t', index=False, header=False)


def parse_args():
    """Parses arguments for pyim-cis."""

    parser = argparse.ArgumentParser(prog='pyim-bed')

    parser.add_argument('--insertions', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)

    parser.add_argument('--width', default=500, type=int)

    parser.add_argument('--drop', nargs='+', default=None)

    return parser.parse_args()


if __name__ == '__main__':
    main()
