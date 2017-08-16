"""Script for the pyim-align command.

The align command is responsible for extracting genomic reads from the
sequencing data, aligning these reads to the reference genome and extracting
insertion sites from these alignments. The command provides access to several
distinct pipelines, which perform these tasks for different types
of sequencing data.
"""

import argparse
import logging

from pyim.align.aligners import AlignerCommand

logging.basicConfig(
    format='[%(asctime)-15s]  %(message)s',
    level=logging.INFO,
    datefmt='%Y-%m-%d %H:%M:%S')


def main():
    """Main function for pyim-align."""

    args = parse_args()
    args.command.run(args)


def parse_args():
    """Parses arguments for pyim-align."""

    # Setup main parser.
    parser = argparse.ArgumentParser(prog='pyim-align')
    subparsers = parser.add_subparsers(dest='aligner')
    subparsers.required = True

    # Register pipelines.
    commands = AlignerCommand.available_commands()

    for name, command in commands.items():
        cmd_parser = subparsers.add_parser(name)
        command.configure(cmd_parser)
        cmd_parser.set_defaults(command=command)

    return parser.parse_args()


if __name__ == '__main__':
    main()
