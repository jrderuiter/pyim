from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
from builtins import (ascii, bytes, chr, dict, filter, hex, input,
                      int, map, next, oct, open, pow, range, round,
                      str, super, zip)
from future.utils import native_str

import logging
from argparse import ArgumentParser
from pathlib import Path

import pandas as pd

from ._logging import print_header, print_footer


def setup_parser():
    parser = ArgumentParser(prog='pyim-merge')

    parser.add_argument('insertions', nargs='+', type=Path)
    parser.add_argument('output', type=Path)

    parser.add_argument('--names', nargs='+', default=None)
    parser.add_argument('--samples', nargs='+', default=None)
    # parser.add_argument('--complement', default=False, action='store_true')

    return parser


def main():
    parser = setup_parser()
    args = parser.parse_args()

    # Get logger and print header.
    logger = logging.getLogger()
    print_header(logger, command='merge')

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
        if name != '':
            frame['id'] = ['{}.{}'.format(name, id_)
                           for id_ in frame['id']]
        ins_frames.append(frame)

    # Merge frames.
    merged = pd.concat(ins_frames, ignore_index=True)

    logger.info('Merging insertions for {} datasets, containing {} samples'
                .format(len(args.insertions), merged['sample'].nunique()))

    # Filter samples if needed.
    if args.samples is not None:
        logger.info('Subsetting dataset to {} samples'
                    .format(len(args.samples)))

        merged_samples = set(merged['sample'])
        for sample in args.samples:
            if sample not in merged_samples:
                logging.warning('- Missing insertions for sample {}'
                                .format(sample))

        mask = merged['sample'].isin(set(args.samples))
        merged = merged.ix[mask]

    # Write output.
    logging.info('Writing merged output')
    merged.to_csv(str(args.output), sep=native_str('\t'), index=False)

    print_footer(logger)


if __name__ == '__main__':
    main()
