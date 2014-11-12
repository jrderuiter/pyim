import argparse

import pandas

from pyim.cis.cimpl import annotate_with_cis


def cis_main(options):
    insertions = pandas.read_csv(options.input, sep='\t')
    cis_insertions = annotate_with_cis(insertions, genome=options.genome,
                                       scales=options.scales, n_iter=options.n_iter,
                                       alpha=options.alpha)
    cis_insertions.to_csv(options.output, sep='\t', index=False, header=True)


def _parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument('-i', '--input', required=True)
    parser.add_argument('-o', '--output', required=True)

    parser.add_argument('-g', '--genome', default='mm10')
    parser.add_argument('-s', '--scales', default=30000, type=int, nargs='+')
    parser.add_argument('-n', '--n_iter', default=1000, type=int)
    parser.add_argument('-a', '--alpha', default=0.05, type=float)

    return parser.parse_args()


if __name__ == '__main__':
    cis_main(_parse_args())
