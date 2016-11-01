from argparse import ArgumentParser
from pathlib import Path

import pandas as pd

from pyim.model import Insertion
from pyim.util import remove_prefix


def setup_parser():
    parser = ArgumentParser(prog='pyim-split')

    parser.add_argument('--insertions', type=Path, required=True)
    parser.add_argument('--output_dir', type=Path, required=True)

    parser.add_argument('--samples', nargs='+', required=False, default=None)
    parser.add_argument('--remove_prefix', default=False, action='store_true')

    return parser


def main():
    parser = setup_parser()
    args = parser.parse_args()

    # Read frame.
    ins_frame = Insertion.from_csv(args.insertions, sep='\t', as_frame=True)

    # Create output directory if it doesn't exist.
    args.output_dir.mkdir(exist_ok=True, parents=True)

    if args.samples is not None:
        # Subset for samples and convert to categorical.
        ins_frame = ins_frame.ix[ins_frame['sample'].isin(args.samples)]
        ins_frame['sample'] = pd.Categorical(
            ins_frame['sample'], categories=args.samples)

    # Split and write individual outputs.
    for sample, sample_frame in ins_frame.groupby('sample'):
        if args.remove_prefix:
            sample_frame['id'] = remove_prefix(
                sample_frame['id'], prefix=sample + '.')

        if len(sample_frame) == 0:
            print('WARNING: no insertions found for sample {}'.format(sample))

        sample_path = args.output_dir / '{}.txt'.format(sample)
        sample_frame.to_csv(str(sample_path), sep='\t', index=False)
