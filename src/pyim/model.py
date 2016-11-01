# pylint: disable=W0622,W0614,W0401
from __future__ import absolute_import, division, print_function
from builtins import *
# pylint: enable=W0622,W0614,W0401

import collections

import pandas as pd
import toolz


class MetadataFrameMixin(object):
    """Mixin class adding namedtuple/frame conversion support."""

    @classmethod
    def _non_metadata_fields(cls):
        fields = list(cls._fields)
        del fields[fields.index('metadata')]
        return fields

    @classmethod
    def to_frame(cls, insertions):
        """Converts list of objects to a dataframe representation."""

        rows = (cls._to_dict(ins) for ins in insertions)

        df = pd.DataFrame.from_records(rows)
        df = cls.format_frame(df)

        return df

    @classmethod
    def format_frame(cls, df):
        """Ensures frame is properly formatted (column order etc.)"""
        cls.check_frame(df)
        return cls._reorder_columns(df, order=cls._non_metadata_fields())

    @classmethod
    def check_frame(cls, df):
        basic_fields = cls._non_metadata_fields()
        missing_columns = set(basic_fields) - set(df.columns)

        if len(missing_columns) > 0:
            raise ValueError('Missing required columns {}',
                             ', '.join(missing_columns))

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
    def from_frame(cls, df):
        """Converts dataframe into a list of objects."""

        cls.check_frame(df)

        basic_fields = cls._non_metadata_fields()
        metadata_fields = list(set(df.columns) - set(basic_fields))

        for row in df.itertuples():
            row_dict = row._asdict()

            metadata = {k: row_dict.pop(k) for k in metadata_fields}
            row_dict.pop('Index', None)

            yield cls(**row_dict, metadata=metadata)


_Insertion = collections.namedtuple(
    'Insertion', ['id', 'chromosome', 'position',
                  'strand', 'metadata'])


class Insertion(MetadataFrameMixin, _Insertion):
    """Model class representing an insertion."""

    __slots__ = ()

    @classmethod
    def format_frame(cls, df):
        df = super().format_frame(df)
        df['position'] = df['position'].astype(int)
        df['strand'] = df['strand'].astype(int)
        return df
