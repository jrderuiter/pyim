"""Script for the pyim-bed command.

Converts an insertion dataframe to the BED file format."""

import argparse
from pathlib import Path

from pyim.model import InsertionSet


def main():
    """Main function for pyim-bed."""

    args = parse_args()

    insertions = InsertionSet.from_csv(args.insertions, sep='\t')
    insertions.to_bed(
        args.output, width=args.width, drop_columns=args.drop_columns)


def parse_args():
    """Parses arguments for pyim-bed."""

    parser = argparse.ArgumentParser(prog='pyim-bed')

    parser.add_argument('--insertions', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)

    parser.add_argument('--width', default=500, type=int)
    parser.add_argument('--drop_columns', nargs='+', default=None)

    return parser.parse_args()


if __name__ == '__main__':
    main()
