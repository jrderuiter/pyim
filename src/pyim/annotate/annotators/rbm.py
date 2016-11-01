from pathlib import Path

from .base import Annotator, CisAnnotator, register_annotator
from .window import WindowAnnotator, Window

# Window format: (us, ua, ds, da)
WINDOW_PRESETS = {
    'SB': (20000, 10000, 25000, 5000),
    'MULV': (20000, 120000, 40000, 5000),
    'MMTV': (20000, 120000, 40000, 5000)
}


class RbmAnnotator(Annotator):
    def __init__(self,
                 reference_gtf,
                 window_sizes=None,
                 preset=None,
                 closest=False,
                 blacklist=None,
                 cis_sites=None):
        super().__init__()

        if window_sizes is None:
            if preset is None:
                raise ValueError('Either window_sizes or '
                                 'preset must be defined')
            else:
                window_sizes = WINDOW_PRESETS[preset]

        windows = self._build_windows(window_sizes)
        self._annotator = WindowAnnotator(
            reference_gtf,
            windows=windows,
            closest=closest,
            blacklist=blacklist)
        self._cis_sites = cis_sites

    @classmethod
    def configure_args(cls, parser):
        super().configure_args(parser)

        # Required arguments.
        parser.add_argument('--reference_gtf', required=True, type=Path)

        # Optional arguments.
        group = parser.add_mutually_exclusive_group(required=True)
        group.add_argument('--preset', choices=WINDOW_PRESETS.keys())
        group.add_argument('--window_sizes', nargs=4, type=int)

        parser.add_argument('--closest', default=False, action='store_true')
        parser.add_argument('--blacklist', nargs='+', default=None)

    @classmethod
    def parse_args(cls, args):
        return {
            'reference_gtf': args.reference_gtf,
            'window_sizes': args.window_sizes,
            'preset': args.preset,
            'closest': args.closest,
            'blacklist': args.blacklist
        }

    def annotate(self, insertions):
        return self._annotator.annotate(insertions)

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

    @property
    def gtf(self):
        return self._annotator.gtf


class RbmCisAnnotator(CisAnnotator, RbmAnnotator):
    pass


register_annotator('rbm', RbmCisAnnotator)

# register_annotator('rbm', RbmAnnotator)
