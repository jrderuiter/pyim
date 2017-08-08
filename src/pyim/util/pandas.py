import itertools
import operator

from natsort import natsorted
import numpy as np
import pandas as pd
import pysam
from intervaltree import IntervalTree


class GenomicDataFrame(pd.DataFrame):
    """DataFrame with fast indexing by genomic position.

    Requires columns 'chromosome', 'start' and 'end' to be present in the
    DataFrame, as these columns are used for indexing.
    """

    _internal_names = pd.DataFrame._internal_names + ['_gi']
    _internal_names_set = set(_internal_names)

    _metadata = ['_chrom_col', '_start_col', '_end_col', '_chrom_lengths']

    def __init__(self,
                 *args,
                 chrom_col='chromosome',
                 start_col='start',
                 end_col='end',
                 chrom_lengths=None,
                 **kwargs):
        super().__init__(*args, **kwargs)

        self._chrom_col = chrom_col
        self._start_col = start_col
        self._end_col = end_col
        self._gi = None

        # TODO: Check lengths for validity.
        self._chrom_lengths = chrom_lengths

    @property
    def gi(self):
        if self._gi is None:
            self._gi = GenomicIndex(
                self,
                chrom_col=self._chrom_col,
                start_col=self._start_col,
                end_col=self._end_col)
        return self._gi

    @property
    def _constructor(self):
        return GenomicDataFrame

    @property
    def chromosomes(self):
        return natsorted(self.gi.chromosome.unique())

    @property
    def chromosome_lengths(self):
        if self._chrom_lengths is None:
            lengths = self.gi.end.groupby(self.gi.chromosome).max()
            self._chrom_lengths = dict(zip(lengths.index, lengths.values))
        return {
            k: v
            for k, v in self._chrom_lengths.items()
            if k in set(self.chromosomes)
        }

    def search(self,
               chromosome,
               begin,
               end,
               strict_left=False,
               strict_right=False):
        """Subsets the DataFrame for rows within given range."""

        return self.gi.search(
            chromosome,
            begin,
            end,
            strict_left=strict_left,
            strict_right=strict_right)

    def with_chromosome_offset(self, suffix='_offset'):
        """Returns copy of DataFrame with offset start/end columns."""

        # Lookup chromosome lengths.
        lengths = pd.Series(self.chromosome_lengths)
        lengths = lengths.loc[self.chromosomes]

        # Lookup offsets.
        offsets = pd.Series(
            np.concatenate([[0], np.cumsum(lengths)]),
            index=self.chromosomes + ['_END_'],
            name='offset')

        df_offsets = offsets.loc[self.gi.chromosome].values

        # Add as new columns.
        offset_start_col = self.gi.start_col + suffix
        offset_end_col = self.gi.end_col + suffix

        new = self.assign(**{
            offset_start_col: self.gi.start + df_offsets,
            offset_end_col: self.gi.end + df_offsets
        })

        new = new.sort_values(by=[offset_start_col, offset_end_col])
        new = new.reindex(new.index).reset_index(drop=True)

        return new, offsets

    @classmethod
    def from_position_df(cls, df, position_col='position', width=0, **kwargs):
        """Builds a GenomicDataFrame from a dataframe with positions."""

        # Note: end is exclusive.
        start_col = kwargs.get('start_col', 'start')
        end_col = kwargs.get('end_col', 'end')

        half_width = width // 2
        df = df.assign(**{
            start_col: df[position_col] - half_width,
            end_col: (df[position_col] - half_width) + 1
        })

        df = df.drop(position_col, axis=1)

        return cls(df, **kwargs)

    @classmethod
    def from_gtf(cls, gtf_path, feature=None):
        """Build a GenomicDataFrame from a GTF file."""

        def _record_to_row(record):
            row = {
                'contig': record.contig,
                'source': record.source,
                'feature': record.feature,
                'start': int(record.start),
                'end': int(record.end),
                'score': record.score,
                'strand': record.strand,
                'frame': record.frame
            }
            row.update(dict(record))
            return row

        # Parse records into rows.
        gtf_file = pysam.TabixFile(str(gtf_path), parser=pysam.asGTF())
        records = (rec for rec in gtf_file.fetch())

        # Filter for specific feature type.
        # TODO: Generalize?
        if feature is not None:
            records = (rec for rec in records if rec.feature == feature)

        # Build dataframe.
        rows = (_record_to_row(rec) for rec in records)
        data = cls(pd.DataFrame.from_records(rows), chrom_col='contig')

        # Reorder columns to correspond with GTF format.
        columns = ('contig', 'source', 'feature', 'start', 'end', 'score',
                   'strand', 'frame')
        data = reorder_columns(data, order=columns)

        return data


class GenomicIndex(object):
    def __init__(self,
                 df,
                 chrom_col='chromosome',
                 start_col='start',
                 end_col='end'):
        self._df = df

        self._chrom_col = chrom_col
        self._start_col = start_col
        self._end_col = end_col

        self._trees = None

    @property
    def df(self):
        return self._df

    @property
    def chromosome(self):
        return self._df[self._chrom_col]

    @property
    def start(self):
        return self._df[self._start_col]

    @property
    def end(self):
        return self._df[self._end_col]

    @property
    def chromosome_col(self):
        return self._chrom_col

    @property
    def start_col(self):
        return self._start_col

    @property
    def end_col(self):
        return self._end_col

    @property
    def trees(self):
        """Returns trees used for indexing the DataFrame."""

        if self._trees is None:
            self._trees = self._build_trees()

        return self._trees

    def _build_trees(self):
        # Subset frame to positions and rename columns to defaults.
        position_df = self._df[
            [self._chrom_col, self._start_col, self._end_col]]
        position_df.columns = ['chromosome', 'start', 'end']

        # Add index and sort by chromosome (for grouping).
        position_df = position_df.assign(index=np.arange(len(self._df)))
        position_df = position_df.sort_values(by='chromosome')

        # Convert to tuples.
        tuples = ((tup.chromosome, tup.start, tup.end, tup.index)
                  for tup in position_df.itertuples())

        return GenomicIntervalTree.from_tuples(tuples)

    def search(self,
               chromosome,
               begin,
               end,
               strict_left=False,
               strict_right=False):
        """Subsets the DataFrame for rows within given range."""

        overlap = self.trees.search(
            chromosome,
            begin,
            end,
            strict_left=strict_left,
            strict_right=strict_right)

        indices = [interval[2] for interval in overlap]

        return self._df.iloc[indices].sort_index()


class GenomicIntervalTree(object):
    """Datastructure for efficiently accessing genomic objects by position."""

    def __init__(self, trees):
        # type: (Dict[str, IntervalTree]) -> None
        self._trees = trees

    def __getitem__(self, i):
        # type: (str) -> IntervalTree
        """Returns tree with given chromosome name."""
        return self._trees[i]

    @classmethod
    def from_tuples(cls, tuples):
        """Builds an instance from tuples.

        Assumes tuples are sorted by chromosome.
        """

        # Group by chromosome.
        groups = itertools.groupby(tuples, key=operator.itemgetter(0))

        # Build trees.
        trees = {}
        for chrom, group in groups:
            trees[chrom] = IntervalTree.from_tuples(
                (start, end, obj) for _, start, end, obj in group)

        return cls(trees)

    def search(self,
               chromosome,
               begin,
               end=None,
               strict_left=False,
               strict_right=False):
        # type: (str, int, int) -> Iterable[object]
        """Searches the tree for objects within given range."""

        overlap = self._trees[chromosome].search(begin, end)

        if strict_left:
            overlap = (int_ for int_ in overlap if int_.begin >= begin)

        if strict_right:
            overlap = (int_ for int_ in overlap if int_.end < end)

        return set(overlap)


def reorder_columns(df, order):
    """Reorders dataframe columns, sorting any extra columns alphabetically."""

    extra_cols = set(df.columns) - set(order)
    return df[list(order) + sorted(extra_cols)]
