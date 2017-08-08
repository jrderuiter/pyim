"""Script for the pyim-bed command.

Converts an insertion dataframe to the BED file format."""

import argparse
from collections import OrderedDict
from pathlib import Path

import numpy as np
import pandas as pd

from pyim.model import Insertion

RED = '255,0,0'
BLUE = '0,0,255'
GRAY = '60,60,60'


def main():
    """Main function for pyim-cis."""

    args = parse_args()

    # Read insertions.
    insertion_df = Insertion.from_csv(args.insertions, sep='\t', as_frame=True)

    # Drop any columns if needed.
    if args.drop_columns is not None:
        insertion_df = insertion_df.drop(args.drop_columns, axis=1)
        insertion_df = insertion_df.drop_duplicates()

    # Convert to BED frame.
    start = (insertion_df['position'] - (args.width // 2)).astype(int)
    end = (insertion_df['position'] + (args.width // 2)).astype(int)

    strand = insertion_df['strand'].map({1: '+', -1: '-', np.nan: '.'})
    color = strand.map({'+': BLUE, '-': RED, '.': GRAY})

    bed_frame = pd.DataFrame(
        OrderedDict([
            ('chrom', insertion_df['chromosome']),
            ('chromStart', start),
            ('chromEnd', end),
            ('name', insertion_df['id']),
            ('score', insertion_df['support']),
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
    parser.add_argument('--drop_columns', nargs='+', default=None)

    return parser.parse_args()


if __name__ == '__main__':
    main()
