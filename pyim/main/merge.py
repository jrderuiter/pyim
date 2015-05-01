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
    parser = ArgumentParser(prog='pyim-merge')

    parser.add_argument('insertions', nargs='+', type=Path)
    parser.add_argument('output', type=Path)

    parser.add_argument('--names', nargs='+', required=False, default=None)

    return parser


def main():
    parser = setup_parser()
    args = parser.parse_args()

    # Generate default names if none given.
    if args.names is None:
        names = ['Set{}'.format(i) for i in range(1, len(args.insertions) + 1)]
    else:
        names = args.names

    # Read frames.
    ins_frames, samples = [], set()
    for (ins_path, name) in zip(args.insertions, names):
        frame = pd.read_csv(str(ins_path), sep=native_str('\t'))

        # Check for overlapping samples.
        frame_samples = set(filter(bool, frame['sample']))
        overlap = samples.intersection(frame_samples)

        if len(overlap) > 0:
            raise ValueError('Overlapping samples between frames ({})'
                             .format(', '.join(overlap)))

        samples = samples.union(frame_samples)

        # Augment ids to avoid duplicates in merged frame.
        frame['insertion_id'] = ['{}.{}'.format(name, id_)
                                 for id_ in frame['insertion_id']]

        ins_frames.append(frame)

    # Merge and write output.
    merged = pd.concat(ins_frames, ignore_index=True)
    merged.to_csv(str(args.output), sep=native_str('\t'), index=False)


if __name__ == '__main__':
    main()
