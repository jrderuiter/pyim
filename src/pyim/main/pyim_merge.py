import logging
from argparse import ArgumentParser
from pathlib import Path
from collections import Counter

import pandas as pd

from pyim.model import Insertion
from pyim.util import add_prefix


def setup_parser():
    parser = ArgumentParser(prog='pyim-merge')

    parser.add_argument('--insertions', nargs='+', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--sample_names', nargs='+', default=None)

    return parser


def main():
    parser = setup_parser()
    args = parser.parse_args()

    # Read and merge frames.
    if args.sample_names is None:
        sample_names = [fp.stem for fp in args.insertions]
    else:
        sample_names = args.sample_names

    # Read frames.
    ins_frames = (pd.read_csv(fp, sep='\t') for fp in args.insertions)

    # Merge and write output.
    merged = merge_frames(ins_frames, sample_names)
    merged.to_csv(str(args.output), sep='\t', index=False)


def merge_frames(insertion_frames, sample_names):
    # Check sample names for duplicates.
    duplicate_samples = [s for s, count in Counter(sample_names).items()
                         if count > 1]

    if len(duplicate_samples) > 1:
        raise ValueError('Duplicate sample names given ({})'
                         .format(', '.join(duplicate_samples)))

    # Merge frames.
    frames = []
    for (frame, sample_name) in zip(insertion_frames, sample_names):
        # Check if frame is valid.
        Insertion.check_frame(frame)

        # Augment frame with sample name.
        frame = frame.copy()
        frame['sample'] = sample_name
        frame['id'] = add_prefix(frame['id'], prefix=sample_name + '.')

        frames.append(frame)

    merged = pd.concat(frames, axis=0)
    merged = Insertion.format_frame(merged)

    return merged


if __name__ == '__main__':
    main()
