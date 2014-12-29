import argparse

import pandas

from pyim.cis.cimpl import annotate_with_cis


def cis_main(args):
    insertions = pandas.read_csv(args.input, sep='\t')

    cis_ins, cis = annotate_with_cis(
        insertions, system=args.system, specificity_pattern=args.specificity_pattern,
        genome=args.genome, chromosomes=args.chromosomes, scales=args.scales,
        n_iter=args.n_iter, alpha=args.alpha, threads=args.threads)

    cis_ins.to_csv(args.output_base + '.cis_ins.txt', sep='\t', index=False, header=True)
    cis.to_csv(args.output_base + '.cis.txt', sep='\t', index=False, header=True)


def _parse_args():
    # TODO: allow self-supplied reference?

    parser = argparse.ArgumentParser()

    parser.add_argument('input')
    parser.add_argument('output_base')

    parser.add_argument('--genome', default='mm10', choices=['mm10'])
    parser.add_argument('--scales', default=30000, type=int, nargs='+')
    parser.add_argument('--n-iter', default=1000, type=int)
    parser.add_argument('--alpha', default=0.05, type=float)
    parser.add_argument('--chromosomes', default=None, nargs='+')
    parser.add_argument('--threads', default=1, type=int)

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--system', choices=['SB', 'PB', 'MMTV', 'MuLV'])
    group.add_argument('--specificity-pattern')

    return parser.parse_args()


def main():
    cis_main(_parse_args())


if __name__ == '__main__':
    main()
