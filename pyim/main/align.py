from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
from builtins import (ascii, bytes, chr, dict, filter, hex, input,
                      int, map, next, oct, open, pow, range, round,
                      str, super, zip)

import argparse
import logging
import pkg_resources

from pyim.pipelines import shear_splink, shear_splink_sb

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
    shear_splink_sb.register(subparsers)

    # Parse args.
    args = parser.parse_args()

    # Dispatch to pipeline.
    version = pkg_resources.require('pyim')[0].version
    header_str = ' PyIM ({}) '.format(version)
    logger.info('{:-^40}'.format(header_str))

    args.main(args)

    logger.info('{:-^40}'.format(' Done! '))


if __name__ == '__main__':
    main()
