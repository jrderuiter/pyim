"""Script for the pyim-align command.

The align command is responsible for extracting genomic reads from the
sequencing data, aligning these reads to the reference genome and extracting
insertion sites from these alignments. The command provides access to several
distinct pipelines, which perform these tasks for different types
of sequencing data.
"""

import argparse
import logging
from pathlib import Path

from pyim.align.aligners import TagmapAligner

from . import Command

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


class AlignerCommand(Command):
    """Base aligner command."""
    pass


class SingleAlignerCommand(AlignerCommand):
    """Base command class for single-end aligners."""

    def configure(self, parser):
        parser.add_argument('--reads', required=True, type=Path)
        parser.add_argument('--output_prefix', required=True, type=Path)


class PairedAlignerCommand(AlignerCommand):
    """Base command class for paired-end aligners."""

    def configure(self, parser):
        parser.add_argument(
            '--reads', nargs=2, required=True, help='Path to input reads.')

        parser.add_argument(
            '--output_prefix',
            required=True,
            type=Path,
            help='Prefix to use for outputs.')


class TagmapCommand(PairedAlignerCommand):
    """Command for the tagmap aligner."""

    name = 'tagmap'

    def configure(self, parser):
        super().configure(parser)

        parser.add_argument(
            '--transposon', required=True, help='Transposon sequence.')
        parser.add_argument(
            '--bowtie_index',
            type=Path,
            required=True,
            help='Path to bowtie index of the reference genome.')

        parser.add_argument(
            '--min_length',
            type=int,
            default=15,
            help='Minimum length for sequences after trimming.')

        parser.add_argument(
            '--min_support',
            type=int,
            default=2,
            help='Minimum number of reads required to support an insertion.')

        parser.add_argument(
            '--min_mapq',
            type=int,
            default=23,
            help='Minimum mapping quality of alignments.')

        parser.add_argument(
            '--merge_distance',
            type=int,
            default=None,
            help='Distance within which alignments are merged '
            'into a single insertion.')

        parser.add_argument(
            '--local',
            default=False,
            action='store_true',
            help='Perform local alignment.')

        parser.add_argument(
            '--threads',
            default=1,
            type=int,
            help='Number of theads to use in trimming/alignment.')

        return parser

    def run(self, args):
        # Create output directory.
        args.output_prefix.parent.mkdir(exist_ok=True, parents=True)

        # Run aligner.
        bowtie_options = {'--local': args.local, '--threads': args.threads}

        aligner = TagmapAligner(
            transposon_seq=args.transposon,
            bowtie_index_path=args.bowtie_index,
            min_length=args.min_length,
            min_support=args.min_support,
            min_mapq=args.min_mapq,
            merge_distance=args.merge_distance,
            bowtie_options=bowtie_options,
            threads=args.threads)

        insertions = aligner.extract_insertions(
            tuple(args.reads), output_prefix=args.output_prefix, verbose=True)

        # Write insertions.
        ins_path = str(args.output_prefix) + '.txt'
        insertions.to_csv(ins_path, sep='\t', index=False)


if __name__ == '__main__':
    main()
