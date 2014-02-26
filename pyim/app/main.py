

import argparse
from pyim.app.main_cimpl import setup_parser as setup_cimpl_parser
from pyim.app.main_gff import setup_parser as setup_gff_parser
from pyim.app.main_merge import setup_parser as setup_merge_parser


def main():
    args = setup_parser().parse_args()
    print args

def setup_parser():
    parser = argparse.ArgumentParser(prog='pyim')

    subparsers = parser.add_subparsers(dest='command')
    setup_cimpl_parser(subparsers)
    setup_gff_parser(subparsers)
    setup_merge_parser(subparsers)

    return parser

if __name__ == '__main__':
    main()