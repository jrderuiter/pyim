
import argparse
import pandas as pd


def main():
    args = setup_parser().parse_args()

    ins_frames =  [pd.read_csv(f, sep='\t', index_col=0) for f in args.input]
    ins_merged = pd.concat(ins_frames, ignore_index=True)

    if args.mask is not None:
        ins_merged = ins_merged[~ins_merged['sample'].isin(args.mask)]

    ins_merged.to_csv(args.output, sep='\t', index=False)


def setup_parser(root=None):
    if root is not None:
        parser = root.add_parser('merge')
    else:
        parser = argparse.ArgumentParser()

    parser.add_argument('-i', '--input', required=True, nargs='+')
    parser.add_argument('-o', '--output', required=True)
    parser.add_argument('-m', '--mask', required=False, nargs='+')

    return parser


if __name__ == '__main__':
    main()