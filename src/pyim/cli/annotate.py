import argparse
from pathlib import Path

from genopandas import GenomicDataFrame

from pyim.annotate.annotators import WindowAnnotator, RbmAnnotator
from pyim.annotate.annotators.rbm import WINDOW_PRESETS as RBM_WINDOW_PRESETS
from pyim.model import InsertionSet

from . import Command


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


class AnnotatorCommand(Command):
    """Base annotator command."""

    def configure(self, parser):
        parser.add_argument('--insertions', type=Path, required=True)
        parser.add_argument('--output', type=Path, required=True)

        parser.add_argument('--gtf', required=True, type=Path)
        parser.add_argument('--closest', default=False, action='store_true')

        parser.add_argument('--blacklist', nargs='+', default=None)
        parser.add_argument('--blackist_col', default='gene_id')

    def run(self, args):
        # Read insertions and genes.
        insertions = InsertionSet.from_csv(args.insertions, sep='\t')

        genes = GenomicDataFrame.from_gtf(
            args.gtf, filter_=lambda rec: rec['feature'] == 'gene')

        # Annotate!
        annotator = self._build_annotator(args)

        annotated = annotator.annotate(
            insertions,
            gene=genes,
            closest=args.closest,
            blacklist=args.blacklist,
            blacklist_col=args.blacklist_col)

        # Write output.
        annotated.to_csv(args.output, sep='\t')

    def _build_annotator(self, args):
        raise NotImplementedError()


class WindowAnnotatorCommand(AnnotatorCommand):
    """WindowAnnotator command."""

    name = 'window'

    def configure(self, parser):
        super().configure(parser)
        parser.add_argument('--window_size', default=20000, type=int)

    def _build_annotator(self, args):
        return WindowAnnotator.from_window_size(window_size=args.window_size)


class RbmAnnotatorCommand(AnnotatorCommand):
    """RbmAnnotator command."""

    name = 'rbm'

    def configure(self, parser):
        super().configure(parser)

        group = parser.add_mutually_exclusive_group(required=True)
        group.add_argument('--preset', choices=RBM_WINDOW_PRESETS.keys())
        group.add_argument('--window_sizes', nargs=4, type=int)

    def _build_annotator(self, args):
        return RbmAnnotator(window_sizes=args.window_sizes, preset=args.preset)


if __name__ == '__main__':
    main()
