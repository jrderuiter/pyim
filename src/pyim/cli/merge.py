import logging
from argparse import ArgumentParser
from collections import Counter
from itertools import chain
from pathlib import Path

import pandas as pd

from pyim.model import Insertion


def main():
    args = parse_args()

    # Read insertions.
    frames = [pd.read_csv(fp, sep='\t') for fp in args.insertions]

    # Check for duplicate samples.
    samples = list(chain.from_iterable(set(df['sample']) for df in frames))
    duplicates = [s for s, count in Counter(samples).items() if count > 1]

    if len(duplicates) > 1:
        raise ValueError('Duplicate sample names in inputs: {}'
                         .format(', '.join(duplicates)))

    # Merge and write output.
    merged = pd.concat(frames, axis=0, ignore_index=True)
    merged.to_csv(str(args.output), sep='\t', index=False)


def parse_args():
    parser = ArgumentParser(prog='pyim-merge')

    parser.add_argument('--insertions', nargs='+', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)

    return parser.parse_args()


if __name__ == '__main__':
    main()
