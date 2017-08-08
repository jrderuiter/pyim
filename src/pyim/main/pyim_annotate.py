import argparse

from natsort import order_by_index, index_natsorted

from pyim.annotate import AnnotatorCommand
from pyim.model import Insertion


def main():
    """Main function for pyim-annotate."""

    args = parse_args()
    args.command.run(args)


def parse_args():
    """Parses arguments for pyim-annotate."""

    # Setup main parser.
    parser = argparse.ArgumentParser(prog='pyim-annotate')
    subparsers = parser.add_subparsers(dest='annotator')
    subparsers.required = True

    # Register pipelines.
    commands = AnnotatorCommand.available_commands()

    for name, command in commands.items():
        cmd_parser = subparsers.add_parser(name)
        command.configure(cmd_parser)
        cmd_parser.set_defaults(command=command)

    return parser.parse_args()


if __name__ == '__main__':
    main()
