from pathlib import Path

from .base import Annotator, AnnotatorCommand, CisAnnotator
from .window import WindowAnnotator, Window
from ..util import select_closest, filter_blacklist

# Window format: (us, ua, ds, da)
WINDOW_PRESETS = {
    'SB': (20000, 10000, 25000, 5000),
    'MULV': (20000, 120000, 40000, 5000),
    'MMTV': (20000, 120000, 40000, 5000)
}


class RbmAnnotator(Annotator):
    def __init__(self,
                 genes,
                 window_sizes=None,
                 preset=None,
                 closest=False,
                 blacklist=None):
        super().__init__()

        if window_sizes is None:
            if preset is None:
                raise ValueError('Either window_sizes or '
                                 'preset must be defined')
            else:
                window_sizes = WINDOW_PRESETS[preset]

        windows = self._build_windows(window_sizes)
        self._annotator = WindowAnnotator(
            windows=windows, genes=genes, closest=closest, blacklist=blacklist)

    def annotate(self, insertions):
        yield from self._annotator.annotate(insertions)

    def _build_windows(self, window_sizes):
        us, ua, ds, da = window_sizes

        windows = [
            Window(0, 1, strand=1, strict_left=False,
                   strict_right=False, name='is'),
            Window(0, 1, strand=-1, strict_left=False,
                   strict_right=False, name='ia'),
            Window(-us, 0, strand=1, strict_left=False,
                   strict_right=True, name='us'),
            Window(-ua, 0, strand=-1, strict_left=False,
                   strict_right=True, name='ua'),
            Window(1, ds, strand=1, strict_left=True,
                   strict_right=False, name='ds'),
            Window(1, da, strand=-1, strict_left=True,
                   strict_right=False, name='da')] # yapf: disable

        return windows


class RbmAnnotatorCommand(AnnotatorCommand):

    name = 'rbm'

    def configure(self, parser):
        super().configure(parser)

        # Required arguments.
        parser.add_argument('--gtf', required=True, type=Path)

        # Optional arguments.
        group = parser.add_mutually_exclusive_group(required=True)
        group.add_argument('--preset', choices=WINDOW_PRESETS.keys())
        group.add_argument('--window_sizes', nargs=4, type=int)

        parser.add_argument('--closest', default=False, action='store_true')
        parser.add_argument('--blacklist', nargs='+', default=None)

        parser.add_argument('--cis_sites', default=None, type=Path)

    def run(self, args):
        # Read insertions and genes.
        insertions = self._read_insertions(args.insertions)
        genes = self._read_genes_from_gtf(args.gtf)

        # Setup annotator.
        if args.cis_sites is not None:
            cis_sites = list(self._read_cis_sites(args.cis_sites))

            sub_annotator = RbmAnnotator(
                genes=genes,
                window_sizes=args.window_sizes,
                preset=args.preset)

            annotator = CisAnnotator(
                annotator=sub_annotator,
                genes=genes,
                cis_sites=cis_sites,
                closest=args.closest,
                blacklist=args.blacklist)
        else:
            annotator = RbmAnnotator(
                genes=genes,
                window_sizes=args.window_sizes,
                preset=args.preset,
                closest=args.closest,
                blacklist=args.blacklist)

        # Annotate insertions and write output.
        annotated = annotator.annotate(insertions)
        self._write_output(args.output, insertions=annotated)
