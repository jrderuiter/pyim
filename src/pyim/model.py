"""Module containing model classes for fusions and insertions."""

from collections import namedtuple

from pyim.util.frozendict import frozendict
import numpy as np
import pandas as pd
import toolz


class MetadataFrameMixin(object):
    """Mixin class adding namedtuple/frame conversion support."""

    _dtypes = {}

    @classmethod
    def _non_metadata_fields(cls):
        fields = list(cls._fields)
        del fields[fields.index('metadata')]
        return fields

    @classmethod
    def to_frame(cls, insertions):
        """Converts list of objects to a dataframe representation."""

        # Check if insertions is empty.
        is_empty, insertions = cls._is_empty(insertions)

        if is_empty:
            df = pd.DataFrame.from_records(
                [], columns=cls._non_metadata_fields())
        else:
            rows = (cls._to_dict(ins) for ins in insertions)
            df = pd.DataFrame.from_records(rows)
            df = cls.format_frame(df)

        return df

    @staticmethod
    def _is_empty(iterable):
        try:
            _, iterable = toolz.peek(iterable)
            empty = False
        except StopIteration:
            empty = True

        return empty, iterable

    @classmethod
    def _to_dict(cls, obj):
        obj_data = obj._asdict()
        metadata = obj_data.pop('metadata')
        return toolz.merge(metadata, obj_data)

    @classmethod
    def _reorder_columns(cls, df, order):
        extra_cols = set(df.columns) - set(order)
        col_order = list(order) + sorted(extra_cols)
        return df[col_order]

    @classmethod
    def check_frame(cls, df):
        missing = set(cls._non_metadata_fields()) - set(df.columns)
        if len(missing) > 0:
            raise ValueError('Missing required columns: {}'
                             .format(', '.join(missing)))

    @classmethod
    def format_frame(cls, df):
        cls.check_frame(df)

        df2 = df.copy()

        for col, dtype in cls._dtypes.items():
            df2[col] = df[col].astype(dtype)

        df2 = cls._reorder_columns(df, order=cls._non_metadata_fields())

        return df2

    @classmethod
    def from_frame(cls, df):
        """Converts dataframe into a list of objects."""

        cls.check_frame(df)

        basic_fields = cls._non_metadata_fields()
        metadata_fields = list(set(df.columns) - set(basic_fields))

        for row in df.itertuples():
            row_dict = row._asdict()

            metadata = {k: row_dict.pop(k) for k in metadata_fields}
            metadata = frozendict(toolz.valfilter(_not_nan, metadata))

            row_dict.pop('Index', None)

            if not set(basic_fields) == set(row_dict.keys()):
                missing_fields = set(basic_fields) - set(row_dict.keys())
                raise ValueError('Missing required fields ({})'
                                 .format(', '.join(missing_fields)))

            yield cls(metadata=metadata, **row_dict)

    @classmethod
    def from_csv(cls, file_path, as_frame=False, **kwargs):
        df = pd.read_csv(file_path, dtype=cls._dtypes, **kwargs)
        cls.check_frame(df)

        if as_frame:
            return df
        else:
            return cls.from_frame(df)

    @classmethod
    def to_csv(cls, file_path, insertions, index=False, **kwargs):
        df = cls.to_frame(insertions)
        df.to_csv(str(file_path), index=index, **kwargs)


_Insertion = namedtuple('Insertion', [
    'id', 'chromosome', 'position', 'strand', 'support', 'sample', 'metadata'
])


class Insertion(MetadataFrameMixin, _Insertion):
    """Model class representing an insertion."""

    __slots__ = ()

    _dtypes = {'chromosome': str}


_CisSite = namedtuple('CisSite',
                      ['id', 'chromosome', 'position', 'strand', 'metadata'])


class CisSite(MetadataFrameMixin, _CisSite):
    """Model class representing an Common Insertion Site (CIS)."""

    __slots__ = ()

    _dtypes = {'chromosome': str}


def _not_nan(value):
    if value is None:
        return False
    elif isinstance(value, str) and value == '':
        return False
    else:
        try:
            return not np.isnan(value)
        except TypeError:
            return True
