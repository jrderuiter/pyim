import argparse
from pathlib import Path

import pandas as pd

from pyim.external.cutadapt import demultiplex_samples


def main():
    """Main function for pyim-demultiplex."""

    args = parse_args()

    # Construct sample mapping if given.
    if args.sample_mapping is not None:
        map_df = pd.read_csv(str(args.sample_mapping), sep='\t')
        sample_mapping = dict(zip(map_df['barcode'], map_df['sample']))
    else:
        sample_mapping = None

    # Perform de-multiplexing.
    demultiplex_samples(
        reads_path=args.reads,
        output_dir=args.output_dir,
        barcode_path=args.barcodes,
        error_rate=args.error_rate,
        sample_mapping=sample_mapping)


def parse_args():
    """Parses arguments for pyim-demultiplex."""

    parser = argparse.ArgumentParser(prog='pyim-demultiplex')

    parser.add_argument('--reads', required=True, type=Path)
    parser.add_argument('--output_dir', required=True, type=Path)
    parser.add_argument('--barcodes', required=True, type=Path)

    parser.add_argument('--sample_mapping', type=Path)
    parser.add_argument('--error_rate', type=float, default=0.0)

    return parser.parse_args()


if __name__ == '__main__':
    main()
