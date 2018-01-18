import argparse
import logging
from pathlib import Path

import pandas as pd

from pyim.model import Insertion

logging.basicConfig(
    format='[%(asctime)-15s]  %(message)s',
    level=logging.INFO,
    datefmt='%Y-%m-%d %H:%M:%S')


def main():
    """Main function for pyim-split."""

    args = parse_args()

    # Read frame.
    insertion_df = Insertion.from_csv(args.insertions, sep='\t', as_frame=True)

    # Create output directory if it doesn't exist.
    args.output_dir.mkdir(exist_ok=True, parents=True)

    if args.samples is not None:
        # Subset for samples and convert to categorical.
        mask = insertion_df['sample'].isin(args.samples)

        insertion_df = insertion_df.loc[mask]
        insertion_df['sample'] = pd.Categorical(
            insertion_df['sample'], categories=args.samples)

    # Split and write individual outputs.
    for sample, grp in insertion_df.groupby('sample'):
        if args.remove_prefix:
            grp['id'] = grp['id'].str.replace(sample + '.', '')

        if len(grp) == 0:
            print('WARNING: no insertions found for sample {}'.format(sample))

        sample_path = args.output_dir / '{}.txt'.format(sample)
        grp.to_csv(str(sample_path), sep='\t', index=False)


def parse_args():
    """Parses arguments for pyim-split."""

    parser = argparse.ArgumentParser(prog='pyim-split')

    parser.add_argument('--insertions', type=Path, required=True)
    parser.add_argument('--output_dir', type=Path, required=True)

    parser.add_argument('--samples', nargs='+', required=False, default=None)
    parser.add_argument('--remove_prefix', default=False, action='store_true')

    return parser.parse_args()


if __name__ == '__main__':
    main()
