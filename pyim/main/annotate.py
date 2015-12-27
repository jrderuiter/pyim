from __future__ import (absolute_import, division,
                        print_function, unicode_literals)

import argparse
import logging

from pyim.annotation import window, rbm, kcrbm
from ._logging import print_header, print_footer


def main():
    logger = logging.getLogger()

    # Setup main parser.
    parser = argparse.ArgumentParser(prog='pyim-annotate')
    subparsers = parser.add_subparsers(dest='annotator')
    subparsers.required = True

    # Register pipelines.
    window.register(subparsers)
    rbm.register(subparsers)
    kcrbm.register(subparsers)

    # Parse args.
    args = parser.parse_args()

    # Dispatch to pipeline.
    cmd_str = '{} {}'.format('annotate', args.annotator)
    print_header(logger, command=cmd_str)
    args.main(args)
    print_footer(logger)


if __name__ == '__main__':
    main()
