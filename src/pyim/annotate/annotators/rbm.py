from .base import Annotator
from .window import WindowAnnotator, Window

# Window format: (us, ua, ds, da)
WINDOW_PRESETS = {
    'SB': (20000, 10000, 25000, 5000),
    'MULV': (20000, 120000, 40000, 5000),
    'MMTV': (20000, 120000, 40000, 5000)
}


class RbmAnnotator(Annotator):
    def __init__(self, window_sizes=None, preset=None):
        super().__init__()

        if window_sizes is None:
            if preset is None:
                raise ValueError('Either window_sizes or '
                                 'preset must be defined')
            else:
                window_sizes = WINDOW_PRESETS[preset]

        windows = self._build_windows(window_sizes)
        self._annotator = WindowAnnotator(windows=windows)

    def _build_windows(self, window_sizes):
        up_sense, up_anti, down_sense, down_anti = window_sizes

        windows = [
            Window(0, 1, strand=1, strict_left=False,
                   strict_right=False, name='is'),
            Window(0, 1, strand=-1, strict_left=False,
                   strict_right=False, name='ia'),
            Window(-up_sense, 0, strand=1, strict_left=False,
                   strict_right=True, name='us'),
            Window(-up_anti, 0, strand=-1, strict_left=False,
                   strict_right=True, name='ua'),
            Window(1, down_sense, strand=1, strict_left=True,
                   strict_right=False, name='ds'),
            Window(1, down_anti, strand=-1, strict_left=True,
                   strict_right=False, name='da')] # yapf: disable

        return windows

    def match(self, insertions, genes):
        """"Matches insertion to predicted target genes."""
        return self._annotator.match(insertions, genes=genes)
