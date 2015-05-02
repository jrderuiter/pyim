from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
from builtins import (ascii, bytes, chr, dict, filter, hex, input,
                      int, map, next, oct, open, pow, range, round,
                      str, super, zip)

from collections import namedtuple

import pytest

import numpy as np
from skbio import DNASequence

from pyim.alignment.vector import ExactAligner
from pyim.pipelines._base import InsertionIdentifier


# -------- Base identifier tests -------

AlignedSegment = namedtuple(
    'AlignedSegment',
    ['query_name', 'reference_id', 'reference_start',
     'reference_end', 'is_reverse'])


@pytest.fixture(scope='module')
def test_data():
    alignments = [
        AlignedSegment('READ1', 0, 20, 70, False),
        AlignedSegment('READ2', 0, 20, 70, False),
        AlignedSegment('READ3', 0, 20, 69, False),
        AlignedSegment('READ4', 0, 20, 68, False),
        AlignedSegment('READ5', 0, 21, 68, False)
    ]

    bc_map = {'READ1': 'BC01',
              'READ2': 'BC01',
              'READ3': 'BC01',
              'READ4': 'BC02',
              'READ5': 'BC02'}

    return alignments, bc_map


@pytest.fixture(scope='module')
def test_data_stranded():
    return [
        AlignedSegment('READ1', 0, 20, 70, False),
        AlignedSegment('READ2', 0, 20, 69, False),
        AlignedSegment('READ3', 0, 20, 69, True),
        AlignedSegment('READ4', 0, 21, 69, True),
        AlignedSegment('READ5', 0, 21, 68, False),
    ]


@pytest.fixture(scope='module')
def test_data_unordered():
    return [
        AlignedSegment('READ1', 0, 20, 70, False),
        AlignedSegment('READ2', 0, 20, 69, False),
        AlignedSegment('READ4', 0, 21, 69, True),
        AlignedSegment('READ5', 0, 21, 68, False),
        AlignedSegment('READ3', 0, 20, 69, True),
    ]

@pytest.fixture(scope='module')
def test_data_first_reverse():
    return [
        AlignedSegment('READ1', 0, 5, 21, True),
        AlignedSegment('READ2', 0, 20, 69, False),
        AlignedSegment('READ3', 0, 22, 69, False),
    ]


# noinspection PyShadowingNames
# noinspection PyMethodMayBeStatic
class TestInsertionIdentifierGroupByPosition(object):

    def test_simple(self, test_data):
        alignments, _ = test_data

        identifier = InsertionIdentifier()
        groups = list(identifier._group_by_position_barcode(alignments))

        # Should identify a total of two groups.
        assert len(groups) == 2

        # Check group 1 for position and membership.
        (pos, strand, bc), alignments = groups[0]
        assert pos == 20
        assert strand == 1
        assert np.isnan(bc)
        assert len(alignments) == 4

    def test_with_barcode(self, test_data):
        alignments, bc_map = test_data

        identifier = InsertionIdentifier()
        groups = list(identifier._group_by_position_barcode(
            alignments, bc_map))

        # Should identify a total of three groups.
        assert len(groups) == 3

        # Check first groups for position and membership.
        for group in groups[0:1]:
            (pos, strand, bc), alignments = group
            assert pos == 20
            assert strand == 1

            assert bc in {'BC01', 'BC02'}

            if bc == 'BC01':
                assert len(alignments) == 3
            else:
                assert len(alignments) == 1

    def test_stranded(self, test_data_stranded):
        alignments = test_data_stranded

        identifier = InsertionIdentifier()
        groups = list(identifier._group_by_position_barcode(alignments))

        # Should identify a total of three groups.
        assert len(groups) == 3

        # Check first group for position and membership.
        (pos, strand, bc), alignments = groups[0]
        assert pos == 20
        assert strand == 1
        assert len(alignments) == 2

        # Check second group for position and membership.
        (pos, strand, bc), alignments = groups[1]
        assert pos == 21
        assert strand == 1
        assert len(alignments) == 1

        # Check third (reverse) group for position and membership.
        (pos, strand, bc), alignments = groups[2]
        assert pos == 69
        assert strand == -1
        assert len(alignments) == 2

    def test_unordered(self, test_data_unordered):
        with pytest.raises(ValueError):
            identifier = InsertionIdentifier()
            list(identifier._group_by_position_barcode(test_data_unordered))

    def test_first_reverse(self, test_data_first_reverse):
        alignments = test_data_first_reverse

        identifier = InsertionIdentifier()
        groups = list(identifier._group_by_position_barcode(alignments))

        # Should have been returned in increasing order, even though
        # the -1 cluster occurs earlier in the alignments.
        positions = [grp[0][0] for grp in groups]
        assert all(np.diff(positions) >= 0)
