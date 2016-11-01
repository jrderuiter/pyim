# pylint: disable=W0622,W0614,W0401
from __future__ import absolute_import, division, print_function
from builtins import *
# pylint: enable=W0622,W0614,W0401

import logging
from argparse import ArgumentParser
from pathlib import Path
from collections import Counter

import pandas as pd

from pyim.model import Insertion
from ._logging import print_header, print_footer


def setup_parser():
    parser = ArgumentParser(prog='pyim-merge')

    parser.add_argument('--insertions', nargs='+', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--sample_names', nargs='+', default=None)

    return parser


def main():
    parser = setup_parser()
    args = parser.parse_args()

    # Get logger and print header.
    logger = logging.getLogger()
    print_header(logger, command='merge')

    # Read and merge frames.
    merge_files(args.insertions, args.output, sample_names=args.sample_names)

    print_footer(logger)


def merge_files(file_paths, output_path, sample_names=None):
    if sample_names is None:
        sample_names = [fp.stem for fp in file_paths]

    ins_frames = (pd.read_csv(fp, sep='\t') for fp in file_paths)

    merged = merge_frames(ins_frames, sample_names)
    merged.to_csv(str(output_path), sep='\t', index=False)


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
        frame['id'] = (sample_name + '.') + frame['id']

        frames.append(frame)

    merged = pd.concat(frames, axis=0)
    merged = Insertion.format_frame(merged)

    return merged


if __name__ == '__main__':
    main()
