from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
from builtins import (ascii, bytes, chr, dict, filter, hex, input,
                      int, map, next, oct, open, pow, range, round,
                      str, super, zip)

from argparse import ArgumentParser
from pathlib import Path

import pysam
import pandas as pd


def setup_parser():
    parser = ArgumentParser(prog='pyim-split')

    parser.add_argument('alignment_bam', type=Path)
    parser.add_argument('read_barcode_map', type=Path)

    parser.add_argument('--output_dir', type=Path, default='.')

    return parser


def main():
    parser = setup_parser()
    args = parser.parse_args()

    # Create output dir.
    if not args.output_dir.exists():
        args.output_dir.mkdir()

    # Read barcodes.
    barcode_map = pd.read_csv(str(args.read_barcode_map), sep='\t')
    barcode_map = dict(zip(barcode_map['read_id'], barcode_map['barcode']))

    # Split reads into separate files.
    with pysam.AlignmentFile(str(args.alignment_bam), 'rb') as in_file:

        out_files = {}
        try:
            # Open output files.
            for sample in set(barcode_map.values()):
                out_name = args.alignment_bam.stem + '.{}.bam'.format(sample)
                out_path = args.output_dir / out_name

                out_files[sample] = pysam.AlignmentFile(
                    str(out_path), 'wb', template=in_file)

            # Write reads to separate files.
            for read in in_file:
                sample = barcode_map[read.query_name]
                out_files[sample].write(read)

        finally:
            for out_path in out_files.values():
                out_path.close()


if __name__ == '__main__':
    main()
