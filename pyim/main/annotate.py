from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
from builtins import (ascii, bytes, chr, dict, filter, hex, input,
                      int, map, next, oct, open, pow, range, round,
                      str, super, zip)
from future.utils import native_str

import argparse

import pandas as pd

from pyim.annotation import RbmAnnotator, WindowAnnotator

ANNOTATORS = {
    'rbm': RbmAnnotator,
    'window': WindowAnnotator
}


def setup_parser():
    parser = argparse.ArgumentParser(prog='pyim-annotate')

    subparsers = parser.add_subparsers(dest='annotator')
    subparsers.required = True

    for name, class_ in ANNOTATORS.items():
        class_.configure_argparser(subparsers, name=name)

    return parser


def main():
    parser = setup_parser()
    args = parser.parse_args()

    # Check if a sub-parser was chosen.
    if args.annotator is None:
        raise ValueError('No annotator was specified as sub-command (choose '
                         'from {})' .format(', '.join(ANNOTATORS.keys())))

    # Parse options and extract main input/output parameters.
    arg_dict = vars(args)

    annotator_name = arg_dict.pop('annotator')
    input_path = arg_dict.pop('input')
    output_path = arg_dict.pop('output')

    # Instantiate chosen annotator and use to annotate input!
    try:
        annotator_class = ANNOTATORS[annotator_name]
    except KeyError:
        raise ValueError('Pipeline \'{}\' does not exist'
                         .format(annotator_name))
    else:
        annotator = annotator_class.from_args(arg_dict)

        in_frame = pd.read_csv(str(input_path), sep=native_str('\t'))
        out_frame = annotator.annotate(in_frame)

        out_frame.to_csv(str(output_path), sep=native_str('\t'), index=False)


if __name__ == '__main__':
    main()
