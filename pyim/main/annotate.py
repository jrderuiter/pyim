from __future__ import (absolute_import, division,
                        print_function, unicode_literals)

import argparse
import logging

from pyim.annotation import window
from ._logging import print_header, print_footer


def main():
    logger = logging.getLogger()

    # Setup main parser.
    parser = argparse.ArgumentParser(prog='pyim-annotate')
    subparsers = parser.add_subparsers(dest='annotator')
    subparsers.required = True

    # Register pipelines.
    window.register(subparsers)

    # Parse args.
    args = parser.parse_args()

    # Dispatch to pipeline.
    print_header(logger, command='annotate')
    args.main(args)
    print_footer(logger)


if __name__ == '__main__':
    main()
