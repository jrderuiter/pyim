from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
from builtins import (ascii, bytes, chr, dict, filter, hex, input,
                      int, map, next, oct, open, pow, range, round,
                      str, super, zip)


import pytest

from tkgeno.io.model import DNASequence
from pyim.alignment.vector import ExactAligner, SswAligner


@pytest.fixture(scope='module')
def test_data():
    reads = [DNASequence('AATGTGTACCAACTGTTG', 'READ1'),  # Query present.
             DNASequence('AATGTGTACCACAGTTTG', 'READ2'),  # Query in reverse.
             DNASequence('AATGTGTACCATAGTTTG', 'READ3')]  # Query missing.

    query = DNASequence('ACTG', 'QUERY')

    return reads, query


# noinspection PyShadowingNames
# noinspection PyMethodMayBeStatic
class TestExactAligner(object):

    def test_simple(self, test_data):
        reads, query = test_data

        aligner = ExactAligner()
        aln = aligner.align(query=query, target=reads[0])

        assert aln.query_id == 'QUERY'
        assert aln.query_start == 0
        assert aln.query_end == 4

        assert aln.target_id == 'READ1'
        assert aln.target_start == 11
        assert aln.target_end == 15
        assert aln.target_strand == 1

        assert aln.type == 'exact'
        assert aln.coverage == 1.0
        assert aln.identity == 1.0

    def test_reverse(self, test_data):
        reads, query = test_data

        aligner = ExactAligner(try_reverse=True)
        aln = aligner.align(query=query, target=reads[1])

        assert aln.target_id == 'READ2'
        assert aln.target_strand == -1
        assert aln.target_start == 11
        assert aln.target_end == 15

        assert aln.query_id == 'QUERY'
        assert aln.query_start == 0
        assert aln.query_end == 4

        assert aln.type == 'exact'
        assert aln.coverage == 1.0
        assert aln.identity == 1.0

    def test_no_reverse(self, test_data):
        reads, query = test_data

        aligner = ExactAligner(try_reverse=False)
        aln = aligner.align(query=query, target=reads[1])

        assert aln is None

    def test_missing(self, test_data):
        reads, query = test_data

        aligner = ExactAligner()
        aln = aligner.align(query=query, target=reads[2])

        assert aln is None


# noinspection PyShadowingNames
# noinspection PyMethodMayBeStatic
class TestSswAligner(object):

    def test_simple(self, test_data):
        reads, query = test_data

        aligner = SswAligner()
        aln = aligner.align(query=query, target=reads[0])

        assert aln.query_id == 'QUERY'
        assert aln.query_start == 0
        assert aln.query_end == 4

        assert aln.target_id == 'READ1'
        assert aln.target_start == 11
        assert aln.target_end == 15
        assert aln.target_strand == 1

        assert aln.type == 'ssw'
        assert aln.coverage == 1.0
        assert aln.identity == 1.0

    def test_reverse(self, test_data):
        reads, query = test_data

        aligner = SswAligner(try_reverse=True)
        aln = aligner.align(query=query, target=reads[1])

        assert aln.target_id == 'READ2'
        assert aln.target_strand == -1
        assert aln.target_start == 11
        assert aln.target_end == 15

        assert aln.query_id == 'QUERY'
        assert aln.query_start == 0
        assert aln.query_end == 4

        assert aln.type == 'ssw'
        assert aln.coverage == 1.0
        assert aln.identity == 1.0