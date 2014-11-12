import argparse
import pandas
from pyim.annotation.kcrbm import KcrbmAnnotator

ANNOTATORS = {
    'kcrbm': KcrbmAnnotator
}


def annotation_main(options):
    ## TODO: generalize arguments to annotator.
    try:
        class_ = ANNOTATORS[options.annotator]
    except KeyError:
        raise KeyError('Unknown annotator "{}"'.format(options.annotator))

    insertions = pandas.read_csv(options.input, sep='\t')
    annotator = class_(genome='mm10', system='SB')

    if options.method == 'gene':
        annotation = annotator.annotate_by_gene(insertions)
    elif options.method == 'transcript':
        annotation = annotator.annotate_by_transcript(insertions)
    else:
        raise ValueError('Unknown annotation method {}'.format(options.method))

    annotation.to_csv(options.output, sep='\t', index=False, header=True)


def _parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument('-i', '--input', required=True)
    parser.add_argument('-o', '--output', required=True)
    parser.add_argument('-m', '--method', required=True)

    parser.add_argument('-a', '--annotator', default='kcrbm')

    return parser.parse_args()


if __name__ == '__main__':
    annotation_main(_parse_args())

