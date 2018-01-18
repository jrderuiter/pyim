import abc
import argparse
import logging
from pathlib import Path

from pyim.cis.callers import CimplCisCaller
from pyim.model import InsertionSet

from . import Command

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


class CisCallerCommand(Command):
    """Base CisCaller command."""

    def configure(self, parser):
        parser.add_argument('--insertions', type=Path, required=True)
        parser.add_argument('--output', type=Path, required=True)
        parser.add_argument(
            '--output_sites', type=Path, required=False, default=None)

    def run(self, args):
        # Read insertions.
        insertions = InsertionSet.from_csv(args.insertions, sep='\t')

        # Call CIS sites.
        caller = self._build_caller(args)
        cis_sites, cis_mapping = caller.call(insertions)

        # Assign CISs to strand based on their insertions.

        # Annotate insertions with CIS IDs.
        insertions_ann = insertions.merge(
            cis_mapping.rename(columns={'insertion_id': 'id'}),
            on='id',
            how='left')

        # Write results.
        if args.output_sites is None:
            cis_path = args.output.with_suffix('.sites.txt')
        else:
            cis_path = args.output_sites

        insertions_ann.to_csv(args.output, sep='\t', index=False)
        cis_sites.to_csv(cis_path, sep='\t', index=False)

    @abc.abstractmethod
    def _build_caller(self, args):
        pass


class CimplCisCallerCommand(CisCallerCommand):

    name = 'cimpl'

    def configure(self, parser):
        super().configure(parser)

        parser.add_argument('--pattern', required=True)

        parser.add_argument('--genome', default='mm10')
        parser.add_argument('--scales', default=(10000, 30000),
                            nargs='+', type=int)  # yapf: disable
        parser.add_argument('--chromosomes', default=None, nargs='+')
        parser.add_argument('--alpha', default=0.05, type=float)
        parser.add_argument('--lhc_method', default='exclude')
        parser.add_argument('--iterations', default=1000, type=int)
        parser.add_argument('--threads', default=1, type=int)
        #parser.add_argument(
        #    '--min_strand_homogeneity', default=0.75, type=float)

    def _build_caller(self, args):
        return CimplCisCaller(
            pattern=args.pattern,
            genome=args.genome,
            chromosomes=args.chromosomes,
            alpha=args.alpha,
            lhc_method=args.lhc_method,
            iterations=args.iterations,
            threads=args.threads)


if __name__ == '__main__':
    main()
