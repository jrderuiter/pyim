from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
from builtins import (ascii, bytes, chr, dict, filter, hex, input,
                      int, map, next, oct, open, pow, range, round,
                      str, super, zip)

import argparse
import logging

from pyim import __version__
from pyim.pipelines import shear_splink

logging.basicConfig(
        format='%(asctime)-15s %(message)s',
        datefmt='[%Y-%m-%d %H:%M:%S]',
        level=logging.INFO)


def main():
    logger = logging.getLogger()

    # Setup main parser.
    parser = argparse.ArgumentParser(prog='pyim-align')

    subparsers = parser.add_subparsers(dest='pipeline')
    subparsers.required = True

    # Register pipelines.
    shear_splink.register(subparsers)

    # Parse args and dispatch.
    header_str = ' PyIM ({}) '.format(__version__)
    logger.info('{:-^40}'.format(header_str))

    args = parser.parse_args()
    args.main(args)

    logger.info('{:-^40}'.format(' Done! '))


if __name__ == '__main__':
    main()
