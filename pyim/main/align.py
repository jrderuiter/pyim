from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
from builtins import (ascii, bytes, chr, dict, filter, hex, input,
                      int, map, next, oct, open, pow, range, round,
                      str, super, zip)

import argparse

from pyim.pipelines.lam_pcr import LamPcrPipeline

PIPELINES = {
    'lam_pcr': LamPcrPipeline
}


def setup_parser():
    parser = argparse.ArgumentParser(prog='pyim-align')

    subparsers = parser.add_subparsers(dest='pipeline')
    subparsers.required = True

    for name, class_ in PIPELINES.items():
        class_.configure_argparser(subparsers, name=name)

    return parser


def main():
    parser = setup_parser()
    args = parser.parse_args()

    # Check if a sub-parser was chosen.
    if args.pipeline is None:
        raise ValueError('No pipeline was specified as sub-command (choose '
                         'from {})' .format(', '.join(PIPELINES.keys())))

    # Parse options and extract main input/output parameters.
    arg_dict = vars(args)

    pipeline_name = arg_dict.pop('pipeline')
    input_path = arg_dict.pop('input')
    output_path = arg_dict.pop('output')

    # Instantiate chosen pipeline and run!
    try:
        pipeline_class = PIPELINES[pipeline_name]
    except KeyError:
        raise ValueError('Pipeline \'{}\' does not exist'.format(pipeline_name))
    else:
        pipeline = pipeline_class.from_args(**arg_dict)
        pipeline.run(input_path, output_path)


if __name__ == '__main__':
    main()
