from __future__ import (absolute_import, division,
                        print_function, unicode_literals)

import argparse
import logging

from pyim.alignment.pipelines import shear_splink, shear_splink_sb
from ._logging import print_header, print_footer


def main():
    logger = logging.getLogger()

    # Setup main parser.
    parser = argparse.ArgumentParser(prog='pyim-align')
    subparsers = parser.add_subparsers(dest='pipeline')
    subparsers.required = True

    # Register pipelines.
    shear_splink.register(subparsers)
    shear_splink_sb.register(subparsers)

    # Parse args.
    args = parser.parse_args()

    # Dispatch to pipeline.
    print_header(logger, command='align')
    args.main(args)
    print_footer(logger)


if __name__ == '__main__':
    main()
