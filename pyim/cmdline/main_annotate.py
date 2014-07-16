
import argparse
import pandas

from pyim.insertions.annotation.kcrbm import kcrbm_assign_genes


def main():
    args = _parse_args()

    insertions = pandas.read_csv(args.insertion_file, sep='\t')
    annotation = kcrbm_assign_genes(insertions)

    annotation.to_csv(args.output, sep='\t', index=False)


def _parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument('-i', '--input', dest='insertion_file', required=True)
    parser.add_argument('-o', '--output', dest='output', required=True)

    return parser.parse_args()


if __name__ == '__main__':
    main()
