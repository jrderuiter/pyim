import argparse

from natsort import index_natsorted

from pyim_common.io import insertion

from pyim.annotation.gff import GffAnnotator
from pyim.annotation.kcrbm import KcrbmAnnotator

ANNOTATORS = {
    'gtf': GffAnnotator,
    'kcrbm': KcrbmAnnotator
}


def add_subparser(subparsers, annotator_class, name):
    parser = subparsers.add_parser(name=name)

    # Add default positional.
    parser.add_argument('input')
    parser.add_argument('output')

    # Add class specific options.
    annotator_class.configure_argparser(parser)

    # Add default options.
    parser.add_argument('--method', default='gene', choices=['gene', 'transcript'])

    return parser


def main():

    # Setup main argument parser and annotator specific sub-parsers.
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(help='sub-command help', dest='annotator')

    for name, class_ in ANNOTATORS.items():
        add_subparser(subparsers, class_, name)

    args = parser.parse_args()

    # Check if a sub-parser was chosen.
    if args.annotator is None:
        raise ValueError('No annotator was specified as sub command '
                         '(choose from {})' .format(', '.join(ANNOTATORS.keys())))

    # Parse options and extract main input/output/method parameters.
    options = vars(args)

    input_path = options.pop('input')
    output_path = options.pop('output')
    method = options.pop('method')

    # Instantiate chosen pipeline and run!
    annotator_name = options.pop('annotator')
    try:
        annotator_class = ANNOTATORS[annotator_name]
    except KeyError:
        raise ValueError('Annotator \'{}\' does not exist'.format(annotator_name))
    else:
        annotator = annotator_class(**options)

        # Load insertions from input file.
        ins_frame = insertion.read_frame(input_path)

        # Annotate insertions using the requested method.
        if method == 'gene':
            annotation = annotator.annotate_by_gene(ins_frame)
        elif method == 'transcript':
            annotation = annotator.annotate_by_transcript(ins_frame)
        else:
            raise ValueError('Unknown method {}'.format(method))

        # Write output annotation.
        annotation = annotation.iloc[index_natsorted(annotation['name'])]
        annotation.to_csv(output_path, sep='\t', index=False)


if __name__ == '__main__':
    main()
