from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
from builtins import (ascii, bytes, chr, dict, filter, hex, input,
                      int, map, next, oct, open, pow, range, round,
                      str, super, zip)

import argparse

# from pyim.pipelines.lam_pcr import LamPcrPipeline
from pyim.pipelines import shear_splink


def main():
    # Setup main parser.
    parser = argparse.ArgumentParser(prog='pyim-align')

    subparsers = parser.add_subparsers(dest='pipeline')
    subparsers.required = True

    # Register pipelines.
    shear_splink.register(subparsers)

    # Parse args and dispatch.
    args = parser.parse_args()
    args.main(args)


if __name__ == '__main__':
    main()
