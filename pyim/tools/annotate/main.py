import sys
import argparse

import pandas as pd

from pyim.tools.annotate.rbm import RbmAnnotator

ANNOTATORS = {
    'rbm': RbmAnnotator
}


def main():
    # Setup main argument parser and annotator specific sub-parsers.
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(help='sub-command help')

    for name, class_ in ANNOTATORS.items():
        class_.register_parser(subparsers, name)

    args = parser.parse_args()

    if 'annotator' not in args:
        parser.print_help()
        sys.exit(2)
    else:
        # Extract input/output parameters.
        options = vars(args)

        input_path = options.pop('input')
        output_path = options.pop('output')

        # Construct annotator.
        class_ = options.pop('annotator')
        annotator = class_(**options)

        # Load input file.
        frame = pd.read_csv(input_path, sep='\t',
                            dtype={'seqname': str, 'location': int, 'strand': int})

        # Do annotation and write outputs!
        result = annotator.annotate(frame)
        result.to_csv(output_path, sep='\t', index=False)


if __name__ == '__main__':
    main()
