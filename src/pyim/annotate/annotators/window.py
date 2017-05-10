import collections
import itertools
from pathlib import Path

from pyim.util.frozendict import frozendict
from tqdm import tqdm
import toolz

from pyim.util.tabix import GtfFile, GtfFrame

from .base import Annotator, CisAnnotator, register_annotator
from ..filter_ import select_closest
from ..metadata import add_metadata
from ..util import build_interval_trees, numeric_strand


class WindowAnnotator(Annotator):
    def __init__(self,
                 reference_gtf,
                 windows,
                 closest=False,
                 blacklist=None,
                 verbose=True):
        super().__init__()

        self._windows = windows
        self._gtf = GtfFile(reference_gtf)
        self._closest = closest
        self._blacklist = blacklist
        self._verbose = verbose

        self._gtf_frame = None
        self._gtf_trees = None

    @classmethod
    def configure_args(cls, parser):
        super().configure_args(parser)

        # Required arguments.
        parser.add_argument('--reference_gtf', required=True, type=Path)

        # Optional arguments.
        parser.add_argument('--window_size', default=20000, type=int)
        parser.add_argument('--closest', default=False, action='store_true')
        parser.add_argument('--blacklist', nargs='+', default=None)

    @classmethod
    def parse_args(cls, args):
        window_size = args.window_size // 2
        windows = [Window(
            -window_size,
            window_size,
            strand=None,
            name=None,
            strict_left=False,
            strict_right=False)]

        return {
            'reference_gtf': args.reference_gtf,
            'windows': windows,
            'closest': args.closest,
            'blacklist': args.blacklist
        }

    def annotate(self, insertions):
        # Annotate insertions.
        if self._verbose:
            insertions = tqdm(list(insertions), ncols=80)

        annotated = itertools.chain.from_iterable(
            (self._annotate_insertion(ins) for ins in insertions))

        # Add metadata.
        annotated = add_metadata(annotated, reference_gtf=self.gtf)

        # Filter to closest if needed.
        if self._closest:
            annotated = list(select_closest(annotated))

        return annotated

    def _annotate_insertion(self, insertion):
        trees = self._trees

        # Identify overlapping features.
        hits = set()
        for window in self._windows:
            applied_window = window.apply(insertion.chromosome,
                                          insertion.position, insertion.strand)

            hits |= {(feature['gene_id'], feature['gene_name'], window.name)
                     for feature in applied_window.get_overlap(trees)}

        # Filter for blacklist.
        if self._blacklist is not None:
            hits = {hit for hit in hits if hit[1] not in self._blacklist}

        if len(hits) > 0:
            # Annotate insertion with overlapping genes.
            for gene_id, gene_name, window_name in hits:
                metadata = {'gene_id': gene_id, 'gene_name': gene_name}

                if window_name is not None:
                    metadata['window'] = window_name

                metadata = toolz.merge(insertion.metadata, metadata)
                yield insertion._replace(metadata=frozendict(metadata))
        else:
            # In case of no overlap, return original insertion.
            yield insertion

    @property
    def gtf(self):
        if self._gtf_frame is None:
            if isinstance(self._gtf, GtfFrame):
                self._gtf_frame = self._gtf
            else:
                self._gtf_frame = GtfFrame.from_records(
                    self._gtf.fetch(filters={'feature': 'gene'}))
        return self._gtf_frame

    @property
    def _trees(self):
        if self._gtf_trees is None:
            self._gtf_trees = build_interval_trees(self.gtf)
        return self._gtf_trees


class WindowCisAnnotator(CisAnnotator, WindowAnnotator):
    pass


register_annotator('window', WindowCisAnnotator)

_Window = collections.namedtuple('Window', ['start', 'end', 'strand', 'name',
                                            'strict_left', 'strict_right'])


class Window(_Window):
    __slots__ = ()

    def apply(self, chromosome, position, strand):
        # Determine start/end position.
        if strand == 1:
            start = position + self.start
            end = position + self.end

            strict_left = self.strict_left
            strict_right = self.strict_right
        elif strand == -1:
            start = position - self.end
            end = position - self.start

            strict_right = self.strict_left
            strict_left = self.strict_right
        else:
            raise ValueError('Unknown value for strand ({})'.format(strand))

        # Determine new strand.
        if self.strand is not None:
            new_strand = self.strand * strand
        else:
            new_strand = None

        return AppliedWindow(chromosome, start, end, new_strand, self.name,
                             strict_left, strict_right)


_AppliedWindow = collections.namedtuple(
    'AppliedWindow', ['chromosome', 'start', 'end', 'strand', 'name',
                      'strict_left', 'strict_right'])


class AppliedWindow(_AppliedWindow):
    __slots__ = ()

    def get_overlap(self, interval_trees):
        # Find overlapping features (end-inclusive).
        try:
            tree = interval_trees[self.chromosome]
            overlap = tree[self.start:self.end + 1]
        except KeyError:
            overlap = []

        # Extract features.
        features = (interval[2] for interval in overlap)

        # Filter inclusive/exclusive if needed.
        if self.strict_left:
            features = (f for f in features if f['start'] > self.start)

        if self.strict_right:
            features = (f for f in features if f['end'] < self.end)

        # Filter for strand if needed.
        if self.strand is not None:
            features = (f for f in features
                        if numeric_strand(f['strand']) == self.strand)

        return features
