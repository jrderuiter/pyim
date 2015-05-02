from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
from builtins import (ascii, bytes, chr, dict, filter, hex, input,
                      int, map, next, oct, open, pow, range, round,
                      str, super, zip)

import pytest

import pandas as pd

from pyim.pipelines._base import genomic_distance
from pyim.cluster import cluster_frame, cluster_frame_merged


@pytest.fixture(scope='module')
def test_data():
    frame = pd.DataFrame(
        [('1', 10, 1, 10), ('1', 5, 1, 12),
         ('1', 40, 1, 20), ('2', 50, 1, 10)],
        columns=['seqname', 'location', 'strand', 'value'])

    def merge_func(grp):
        ref = grp.iloc[0]
        return pd.Series({
            'seqname': ref.seqname,
            'location': int(grp.location.mean()),
            'strand': ref.strand,
            'value': grp.value.sum()
        }, index=grp.columns)

    dist_func = genomic_distance

    return frame, dist_func, merge_func


# noinspection PyShadowingNames
# noinspection PyMethodMayBeStatic
class TestClusterFrameMerged(object):

    def test_simple(self, test_data):
        frame, dist_func, merge_func = test_data

        res = cluster_frame_merged(frame, dist_func, merge_func, t=10)

        assert isinstance(res, pd.DataFrame)
        assert len(res) == 2

        # Should have grouped the first two entries.
        assert tuple(res.iloc[0]) == ('1', 7, 1, 22)

        # Due to omission of groupby, merge should naively have
        # merged the last two entries into a single one, using
        # the non-summarized entries (seqname, strand)
        # from the first row.
        assert tuple(res.iloc[1]) == ('1', 45, 1, 30)

    def test_grouping(self, test_data):
        frame, dist_func, merge_func = test_data

        res = cluster_frame_merged(frame, dist_func, merge_func,
                                   groupby=['seqname'], t=10)

        assert isinstance(res, pd.DataFrame)
        assert len(res) == 3

        # Should have grouped the first two entries.
        assert tuple(res.iloc[0]) == ('1', 7, 1, 22)

        # Should not have grouped the last entries,
        # due to the difference in seqname.
        assert tuple(res.iloc[1]) == ('1', 40, 1, 20)
        assert tuple(res.iloc[2]) == ('2', 50, 1, 10)

    def test_single(self, test_data):
        frame, dist_func, merge_func = test_data

        res = cluster_frame_merged(frame.iloc[0:1], dist_func, merge_func,
                                   groupby=['seqname'], t=10)

        assert isinstance(res, pd.DataFrame)
        assert len(res) == 1
        assert tuple(res.iloc[0]) == ('1', 10, 1, 10)

    def test_empty(self, test_data):
        frame, dist_func, merge_func = test_data

        frame = pd.DataFrame.from_records([], columns=frame.columns)
        res = cluster_frame_merged(frame, dist_func, merge_func,
                                   groupby=['seqname'], t=10)

        assert isinstance(res, pd.DataFrame)
        assert len(res) == 0
        assert all(res.columns == frame.columns)

    def test_invalid_linkage(self, test_data):
        frame, dist_func, merge_func = test_data

        with pytest.raises(ValueError):
            cluster_frame_merged(frame, dist_func, merge_func,
                                 linkage='whatever', t=10)
