
import argparse
import pandas as pd
from os import path

from pyim.insertions import write_to_gff

def main():
    args = setup_parser().parse_args()
    ins_frame = pd.read_csv(args.input, sep='\t')
    output = path.splitext(args.input)[0] + '.gff' if args.output is None else args.output
    write_to_gff(ins_frame, output)


def setup_parser(root=None):
    if root is not None:
        parser = root.add_parser('gff')
    else:
        parser = argparse.ArgumentParser()

    parser.add_argument('-i', '--input', required=True)
    parser.add_argument('-o', '--output', required=False, default=None)

    return parser


if __name__ == '__main__':
    main()