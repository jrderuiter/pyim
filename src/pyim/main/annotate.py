from __future__ import absolute_import, division, print_function

#pylint: disable=wildcard-import,unused-wildcard-import,redefined-builtin
from builtins import *
#pylint: enable=wildcard-import,unused-wildcard-import,redefined-builtin

import argparse

import pandas as pd

from pyim.annotation import get_annotators
from pyim.model import Insertion

# pylint: disable=import-error
# from ._logging import print_header, print_footer
# pylint: enable=import-error


def main():
    args = parse_args()

    ins_frame = pd.read_csv(args.input, sep='\t')
    insertions = Insertion.from_frame(ins_frame)

    annotator = args.class_.from_args(args)
    annotated = list(annotator.annotate(insertions))

    annotated_frame = Insertion.to_frame(annotated)
    annotated_frame.to_csv(args.output, sep='\t', index=False)

    # Dispatch to pipeline.
    #cmd_str = '{} {}'.format('annotate', args.annotator)
    #print_header(logger, command=cmd_str)
    #args.main(args)
    #print_footer(logger)


def parse_args():
    # Setup main parser.
    parser = argparse.ArgumentParser(prog='pyim-annotate')
    subparsers = parser.add_subparsers(dest='annotator')
    subparsers.required = True

    # Register pipelines.
    for name, class_ in get_annotators().items():
        annot_parser = subparsers.add_parser(name)
        
        _add_default_arguments(annot_parser)
        class_.setup_args(annot_parser)
        
        annot_parser.set_defaults(class_=class_)
        
    # Actually parse args.
    args = parser.parse_args()

    return args


def _add_default_arguments(parser):
    parser.add_argument('input')
    parser.add_argument('output')


if __name__ == '__main__':
    main()
