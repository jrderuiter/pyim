import argparse
import pandas

from pyim.io import write_insertions_to_gff


def main():
    args = _parse_args()
    insertions = pandas.read_csv(args.insertions_file, sep='\t')
    write_insertions_to_gff(insertions, args.output_file)


def _parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument('-i', '--input_file', dest='insertions_file', required=True)
    parser.add_argument('-o', '--output_file', dest='output_file', required=True)

    return parser.parse_args()


if __name__ == '__main__':
    main()
