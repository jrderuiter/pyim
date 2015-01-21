import argparse

from pyim.alignment.pipelines.lam_pcr import LamPcrPipeline

PIPELINES = {
    'lam_pcr': LamPcrPipeline
}


def add_subparser(subparsers, annotator_class, name):
    parser = subparsers.add_parser(name=name)

    # Add default options.
    parser.add_argument('input')
    parser.add_argument('output')

    # Add class specific options.
    annotator_class.configure_argparser(parser)

    return parser


def main():

    # Setup main argument parser and annotator specific sub-parsers.
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(help='sub-command help', dest='pipeline')

    for name, class_ in PIPELINES.items():
        add_subparser(subparsers, class_, name)

    args = parser.parse_args()

    # Check if a sub-parser was chosen.
    if args.pipeline is None:
        raise ValueError('No pipeline was specified as sub-command '
                         '(choose from {})' .format(', '.join(PIPELINES.keys())))

    # Parse options and extract main input/output parameters.
    options = vars(parser.parse_args())
    input_path = options.pop('input')
    output_path = options.pop('output')

    # Instantiate chosen pipeline and run!
    pipeline_name = options.pop('pipeline')
    try:
        pipeline_class = PIPELINES[pipeline_name]
    except KeyError:
        raise ValueError('Pipeline \'{}\' does not exist'.format(pipeline_name))
    else:
        pipeline = pipeline_class(**options)
        pipeline.run(input_path, output_path)


if __name__ == '__main__':
    main()