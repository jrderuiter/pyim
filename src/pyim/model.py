from collections import namedtuple, OrderedDict

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


class RecordSet(object):
    """Base class that provides functionality for serializing and
       deserializing namedtuple records into a DataFrame format.

    Subclasses should override the ``_tuple_class`` method to
    return the namedtuple class that should be used as a record.
    """

    def __init__(self, values: pd.DataFrame) -> None:
        self._values = self._check_frame(values)

    @classmethod
    def _check_frame(cls, values):
        fields = cls._tuple_fields()

        for field in fields:
            if field not in values.columns:
                raise ValueError('Missing required column {}'.format(field))

        return values.reindex(columns=fields)

    @classmethod
    def _tuple_class(cls):
        """Returns namedtuple class used to instantiate records."""
        raise NotImplementedError()

    @classmethod
    def _tuple_fields(cls):
        """Returns the fields in the named tuple class."""
        return cls._tuple_class()._fields

    @property
    def values(self) -> pd.DataFrame:
        """Internal DataFrame representation of records."""
        return self._values

    def __getitem__(self, item):
        return self._values[item]

    def __setitem__(self, idx, value):
        self._values[idx] = value

    def __len__(self):
        return len(self._values)

    @classmethod
    def from_tuples(cls, tuples):
        """Builds a record set instance from the given tuples."""
        records = (tup._asdict() for tup in tuples)
        return cls(pd.DataFrame.from_records(records))

    def to_tuples(self):
        """Converts the record set into an iterable of tuples."""

        tuple_class = self._tuple_class()

        for row in self._values.itertuples():
            row_dict = row._asdict()
            row_dict.pop('Index', None)

            yield tuple_class(**row_dict)

    @classmethod
    def from_csv(cls, file_path, **kwargs):
        """Reads a record set from a csv file using pandas.read_csv."""
        values = pd.read_csv(file_path, **kwargs)
        return cls(values)

    def to_csv(self, file_path, **kwargs):
        """Writes the record set to a csv file using pandas' to_csv."""
        self._values.to_csv(file_path, **kwargs)

    def groupby(self, by, **kwargs):
        """Groups the set by values of the specified columns."""
        for key, group in self._values.groupby(by, **kwargs):
            yield key, self.__class__(group)

    def query(self, expr, **kwargs):
        """Queries the columns of the set with a boolean expression."""
        return self.__class__(self._values.query(expr, **kwargs))

    @classmethod
    def concat(cls, record_sets):
        """Concatenates multiple records sets into a single set."""
        return cls(pd.concat((rs.values for rs in record_sets), axis=0))


class MetadataRecordSet(RecordSet):
    """Base RecordSet that supports record metadata.

    Extension of the RecordSet class, which assumes that records contain
    a dict 'metadata' field which contains variable metadata. The
    MetadataRecordSet class ensures that this data is expanded from the
    original record when converted to the set's DataFrame format, and
    converted back again when transforming back to tuples.
    """

    METADATA_FIELD = 'metadata'

    @property
    def metadata_columns(self):
        """Available metadata columns."""
        return set(self._values.columns) - set(self._tuple_fields())

    @classmethod
    def _check_frame(cls, values):
        fields = [
            field for field in cls._tuple_fields()
            if field != cls.METADATA_FIELD
        ]

        for field in fields:
            if field not in values.columns:
                raise ValueError('Missing required column {}'.format(field))

        extra_cols = set(values.columns) - set(fields)
        col_order = list(fields) + sorted(extra_cols)

        return values.reindex(columns=col_order)

    @classmethod
    def from_tuples(cls, tuples):
        """Builds a record set instance from the given tuples."""

        metadata_field = cls.METADATA_FIELD

        def _to_record(tup):
            record = tup._asdict()
            record.update(record.pop(metadata_field))
            return record

        records = (_to_record(tup) for tup in tuples)
        return cls(pd.DataFrame.from_records(records))

    def to_tuples(self):
        """Converts the record set into an iterable of tuples."""

        tuple_class = self._tuple_class()

        # Determine metadata fields.
        tuple_fields = set(self._tuple_fields()) - {self.METADATA_FIELD}
        metadata_fields = set(self._values.columns) - tuple_fields

        for row in self._values.itertuples():
            row_dict = row._asdict()
            row_dict.pop('Index', None)

            # Extract metadata into separate dict.
            metadata = {k: row_dict.pop(k) for k in metadata_fields}

            yield tuple_class(metadata=metadata, **row_dict)


Insertion = namedtuple('Insertion', [
    'id', 'chromosome', 'position', 'strand', 'support', 'sample', 'metadata'
])


class InsertionSet(MetadataRecordSet):
    """Class that represents an insertion dataset."""

    def __init__(self, values, gene_col='gene_name'):
        super().__init__(values)
        self._gene_col = gene_col

    @classmethod
    def _tuple_class(cls):
        return Insertion

    @property
    def samples(self):
        """List of samples in the InsertionSet."""
        return set(self._values['sample'].unique())

    @property
    def genes(self):
        """List of samples in the InsertionSet."""
        return set(self._values[self._gene_col].dropna().unique())

    def subset(self, genes=None, samples=None):
        """Subsets the set to the given genes/samples."""

        values = self._values

        if genes is not None:
            values = values.loc[values[self._gene_col].isin(set(genes))]

        if samples is not None:
            values = values.loc[values['sample'].isin(set(samples))]

        return self.__class__(values.copy())

    def rename(self, name_map, drop=False):
        """Renames samples using the given mapping."""

        if drop:
            mask = self._values['sample'].isin(name_map)
            values = self._values.loc[mask].copy()
        else:
            values = self._values.copy()

        values['sample'] = values['sample'].map(name_map)

        return self.__class__(values)

    def with_relative_support(self, col_name='support_relative'):
        """Adds a relative support score (scaled for the max support in the
           corresponding sample) to the record dataframe.
        """

        def _add_relative(grp):
            rel_support = grp['support'] / grp['support'].max()
            return grp.assign(**{col_name: rel_support})

        values = pd.concat(
            (_add_relative(grp) for _, grp in self._values.groupby('sample')),
            axis=0)

        return self.__class__(values)

    def as_matrix(self,
                  index=None,
                  index_order=None,
                  value='support',
                  fill_value=None,
                  aggfunc='max',
                  sort=False):
        """Returns a index-by-sample value matrix of the insertions."""

        if index is None:
            index = self._gene_col

        matrix = pd.pivot_table(
            self._values,
            index=index,
            columns='sample',
            values=value,
            aggfunc=aggfunc,
            fill_value=fill_value)

        if index_order is not None:
            matrix = matrix.reindex(index=index_order)

        if sort:
            matrix = sort_matrix(matrix, index_order=index_order)

        return matrix

    def plot_matrix(self,
                    index=None,
                    index_order=None,
                    value='support',
                    ax=None,
                    sort=True,
                    heatmap_kws=None):
        """Plots a index-by-sample value matrix of insertions as a heatmap."""

        matrix = self.as_matrix(
            index=index, value=value, index_order=index_order, sort=sort)

        cmap = sns.light_palette('#376CA3', as_cmap=True)
        ax = sns.heatmap(
            matrix,
            cmap=cmap,
            ax=ax,
            cbar_kws={'label': value},
            **(heatmap_kws or {}))

        plt.setp(ax.get_xticklabels(), rotation=90)
        plt.setp(ax.get_yticklabels(), rotation=0)

        for loc in ['top', 'right', 'left', 'bottom']:
            ax.spines[loc].set_visible(True)

        return ax

    def plot_matrix_grouped(self,
                            group,
                            group_order=None,
                            index=None,
                            index_order=None,
                            value='support',
                            figsize=None,
                            gridspec_kws=None,
                            heatmap_kws=None,
                            **kwargs):
        """Plots a grouped index-by-sample value matrix of insertions."""

        # Group samples.
        if group_order is None:
            group_order = list(self._values[group].unique())

        grouped = dict(self.groupby(group))

        # Setup axes.
        width_ratios = [len(grouped[grp].samples) for grp in group_order]

        default_gridspec_kws = {'width_ratios': width_ratios, 'wspace': 0.1}
        gridspec_kws = {**default_gridspec_kws, **(gridspec_kws or {})}

        # Determine vmin/vmax limits.
        vmin, vmax = self._calc_value_limits(
            value, index=index, index_order=index_order)

        default_heatmap_kws = {'vmin': vmin, 'vmax': vmax}
        heatmap_kws = {**default_heatmap_kws, **(heatmap_kws or {})}

        # Draw matrices.
        fig, axes = plt.subplots(
            ncols=len(grouped), figsize=figsize, gridspec_kw=gridspec_kws)

        for (grp_key, ax), last in lookahead(zip(group_order, axes)):
            group_values = grouped[grp_key]

            grp_heatmap_kws = dict(heatmap_kws)
            grp_heatmap_kws['cbar'] = last

            group_values.plot_matrix(
                index=index,
                index_order=index_order,
                value=value,
                ax=ax,
                heatmap_kws=grp_heatmap_kws,
                **kwargs)

            ax.set_title(grp_key)

        # Clean-up axes.
        for ax in axes[1:]:
            ax.set_ylabel('')
            ax.set_yticks([])

        for ax in axes:
            ax.set_xlabel('')

        return fig

    def _calc_value_limits(self, value, index=None, index_order=None):
        """Calculates vmin/vmax for given value and index values."""

        if index is None:
            index = self._gene_col

        if index_order is not None:
            subset = self._values.loc[self._values[index].isin(index_order)]
        else:
            subset = self._values

        vmin = subset[value].min()
        vmax = subset[value].max()

        return vmin, vmax

    def gene_frequency(self):
        """Calculates the frequency (= number of samples) of each gene."""
        frequencies = self._values.groupby(self._gene_col)['sample'].nunique()
        return frequencies.sort_values(ascending=False)

    def gene_sense_fraction(self):
        """Calculates the fraction of sense insertions per gene."""

        # TODO: column check.

        def _sense_frac(x):
            is_sense = x['gene_orientation'] == 'sense'
            return np.sum(is_sense) / len(x)

        return self._values.groupby(self._gene_col).apply(_sense_frac)

    def gene_sense_fraction_weighted(self, weight='support'):
        """Calculates the weighted fraction of sense insertions per gene."""

        def _sense_frac_weighted(x):
            is_sense = x['gene_orientation'] == 'sense'
            return np.sum(is_sense * x[weight]) / np.sum(x[weight])

        return self._values.groupby(self._gene_col).apply(_sense_frac_weighted)

    def sample_frequency(self, collapse_gene=False):
        """Calculates the frequency (= number of insertions) of each sample."""

        if collapse_gene:
            frequency = self._values.groupby('sample')[
                self._gene_col].nunique()
        else:
            frequency = self._values.groupby('sample').size()

        return frequency.sort_values(ascending=False)

    def sample_genes(self, sample):
        mask = self._values['sample'] == sample
        return set(self._values.loc[mask, self._gene_col].dropna())

    def compare_samples(self, sample_a, sample_b, value='support'):
        """Compares gene insertions between two samples."""

        values = self.as_matrix(value=value)[[sample_a, sample_b]]
        return values.dropna(how='all')

    def plot_sample_correlation(self,
                                sample_a,
                                sample_b,
                                value='support',
                                keep_missing=True,
                                ax=None,
                                **kwargs):
        """Plots correlation between two samples by comparing the
           depths of insertions in (matching) genes.
        """

        if ax is None:
            _, ax = plt.subplots()

        plot_data = self.as_matrix(value=value)[[sample_a, sample_b]]

        if keep_missing:
            plot_data = plot_data.fillna(0)
        else:
            plot_data = plot_data.dropna()

        sns.regplot(data=plot_data, x=sample_a, y=sample_b, ax=ax, **kwargs)

        return ax

    def plot_sample_overlap(self, samples, labels=None, ax=None, **kwargs):
        """Plots overlap between two or three samples (in terms of genes)."""

        if ax is None:
            _, ax = plt.subplots()

        if len(samples) == 2:
            from matplotlib_venn import venn2 as venn
        elif len(samples) == 3:
            from matplotlib_venn import venn3 as venn
        else:
            raise ValueError('Only plotting for 2 or 3 samples is supported')

        genes = [self.sample_genes(sample) for sample in samples]

        return venn(genes, set_labels=labels or samples, ax=ax, **kwargs)

    def plot_gene_support(self,
                          order=None,
                          top_n=10,
                          value='support',
                          max_per_sample=False,
                          sample_design=None,
                          ax=None,
                          **kwargs):
        """Plots boxplot of support for genes across samples."""

        if max_per_sample:
            # Take only the strongest insertion per sample (for each gene).
            plot_data = self._values.groupby(
                ['sample', self._gene_col])[value].max().reset_index()
        else:
            plot_data = self._values

        if sample_design is not None:
            sample_design = sample_design.copy()
            sample_design.index.name = 'sample'

            plot_data = pd.merge(
                plot_data,
                sample_design.reset_index(),
                on='sample',
                how='left')

        if order is None:
            # If not order is given, we limit ourselves to the top n genes.
            order = list(self.gene_frequency().head(n=top_n).index)

        ax = sns.boxplot(
            data=plot_data, x=self._gene_col, y=value, order=order, **kwargs)
        plt.setp(ax.get_xticklabels(), rotation=90)

        return ax

    def to_bed(self, file_path, width=500, drop_columns=None):
        """Writes insertions as bed file.

        Useful for external viewers such as IGV.
        """

        values = self._values

        # Remove specified columns and drop resulting duplicates. Prevents
        # writing the same insertion with slightly different annotations
        # multiple times (e.g. same insertion with diff. genes).
        if drop_columns is not None:
            values = values.drop(drop_columns, axis=1)

        values = values.drop_duplicates()

        # Convert to BED frame.
        half_width = width // 2
        start = (values['position'] - half_width).astype(int)
        end = (values['position'] + half_width).astype(int)

        strand = values['strand'].map({1: '+', -1: '-', np.nan: '.'})
        color = strand.map({'+': '0,0,255', '-': '255,0,0', '.': '60,60,60'})

        bed_data = pd.DataFrame(
            OrderedDict([
                ('chrom', values['chromosome']),
                ('chromStart', start),
                ('chromEnd', end),
                ('name', values['id']),
                ('score', values['support']),
                ('strand', strand),
                ('thickStart', start),
                ('thickEnd', end),
                ('itemRgb', color)
            ])
        )  # yapf: disable

        # Write output.
        bed_data.to_csv(file_path, sep='\t', index=False, header=False)


def sort_matrix(matrix, index_order=None):
    """Sorts a 2D matrix, first by its rows and then by its columns.

    Parameters
    ----------
    matrix : pd.DataFrame
        Matrix to sort.
    row_order : List[str]
        Fixed order to keep the rows in. If not given, rows
        are sorted by their number of non-zero values.

    Returns
    -------
    pd.DataFrame
        Sorted matrix.
    """

    if index_order is None:
        freqs = (matrix > 0).sum(axis=1)
        index_order = list(freqs.sort_values(ascending=False).index)

    matrix_sorted = matrix.loc[index_order]
    matrix_sorted = matrix_sorted.T.sort_values(
        by=index_order, ascending=False).T

    return matrix_sorted


def corrcoef(matrix):
    """Calculates the correlation coef for a DataFrame using np.corrcoef."""

    corr = np.corrcoef(matrix)

    if isinstance(matrix, pd.DataFrame):
        corr = pd.DataFrame(corr, index=matrix.index, columns=matrix.index)

    return corr


def lookahead(iterable):
    """Pass through all values from the given iterable, augmented by the
    information if there are more values to come after the current one
    (False), or if it is the last value (True).
    """
    # Get an iterator and pull the first value.
    iterator = iter(iterable)
    last = next(iterator)

    # Run the iterator to exhaustion (starting from the second value).
    for val in iterator:
        # Report the *previous* value (more to come).
        yield last, False
        last = val

    # Report the last value.
    yield last, True
