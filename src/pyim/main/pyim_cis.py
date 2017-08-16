import argparse
import logging

import toolz

from pyim.cis import CisCallerCommand
from pyim.model import Insertion, CisSite
from pyim.vendor.frozendict import frozendict


logging.basicConfig(
    format='[%(asctime)-15s]  %(message)s',
    level=logging.INFO,
    datefmt='%Y-%m-%d %H:%M:%S')


def main():
    """Main function for pyim-cis."""

    args = parse_args()
    args.command.run(args)


def parse_args():
    """Parses arguments for pyim-cis."""

    # Setup main parser.
    parser = argparse.ArgumentParser(prog='pyim-cis')
    subparsers = parser.add_subparsers(dest='caller')
    subparsers.required = True

    # Register pipelines.
    commands = CisCallerCommand.available_commands()

    for name, command in commands.items():
        cmd_parser = subparsers.add_parser(name)
        command.configure(cmd_parser)
        cmd_parser.set_defaults(command=command)

    return parser.parse_args()


if __name__ == '__main__':
    main()
