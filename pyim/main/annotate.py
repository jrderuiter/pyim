from __future__ import absolute_import, division, print_function

#pylint: disable=wildcard-import,unused-wildcard-import,redefined-builtin
from builtins import *
#pylint: enable=wildcard-import,unused-wildcard-import,redefined-builtin

import argparse
import logging

from pyim.annotation.annotator import window, rbm, rbm_cis

# pylint: disable=import-error
from ._logging import print_header, print_footer
# pylint: enable=import-error


def main():
    logger = logging.getLogger()

    # Setup main parser.
    parser = argparse.ArgumentParser(prog='pyim-annotate')
    subparsers = parser.add_subparsers(dest='annotator')
    subparsers.required = True

    # Register pipelines.
    window.register(subparsers)
    rbm.register(subparsers)
    rbm_cis.register(subparsers)
    # kcrbm.register(subparsers)

    # Parse args.
    args = parser.parse_args()

    # Dispatch to pipeline.
    cmd_str = '{} {}'.format('annotate', args.annotator)
    print_header(logger, command=cmd_str)
    args.main(args)
    print_footer(logger)


if __name__ == '__main__':
    main()
