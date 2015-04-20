from collections import namedtuple

import pytest

import numpy as np
from skbio import DNASequence

from pyim.pipelines.base import BasicGenomicExtractor, InsertionIdentifier


# -------- Extractor tests -------

@pytest.fixture(scope='module')
def reads():
    return [
        DNASequence('ACTGAAATGCGTCTGCCCC'),  # BC01, linker
        DNASequence('ACTGAAAAGCGTCTGCCCC'),  # BC02, linker
        DNASequence('ACTGAAATGCGTCTG'),      # BC01, no linker
        DNASequence('ACTGGCGTCTG')       # No bc, no linker
    ]


@pytest.fixture(scope='module')
def transposon_seq():
    return DNASequence('ACTG', id='transposon')


@pytest.fixture(scope='module')
def barcode_seqs():
    return [DNASequence('AAAT', 'BC01'),
            DNASequence('AAAA', 'BC02')]


@pytest.fixture(scope='module')
def linker_seq():
    return DNASequence('CCCC', id='linker')


# noinspection PyShadowingNames
# noinspection PyMethodMayBeStatic
class TestBasicGenomicExtractor(object):

    def test_exact(self, reads, transposon_seq, barcode_seqs, linker_seq):
        extractor = BasicGenomicExtractor(
            transposon_seq, barcode_seqs,
            barcode_map=None, linker_sequence=linker_seq)

        genomic, barcode = extractor.extract_read(reads[0])
        assert genomic.sequence == 'GCGTCTG'
        assert barcode == 'BC01'

        genomic, barcode = extractor.extract_read(reads[1])
        assert genomic.sequence == 'GCGTCTG'
        assert barcode == 'BC02'

    def test_exact_no_linker(self, reads, transposon_seq, barcode_seqs):
        # Instance without linker sequence, should extract.
        extractor = BasicGenomicExtractor(
            transposon_seq, barcode_seqs,
            barcode_map=None, linker_sequence=None)

        genomic, barcode = extractor.extract_read(reads[2])
        assert genomic.sequence == 'GCGTCTG'
        assert barcode == 'BC01'

    def test_exact_no_linker_neg(self, reads, transposon_seq,
                                 barcode_seqs, linker_seq):
        # Instance with linker sequence, should return None.
        extractor = BasicGenomicExtractor(
            transposon_seq, barcode_seqs,
            barcode_map=None, linker_sequence=linker_seq)

        res = extractor.extract_read(reads[2])
        assert res is None

    def test_exact_no_barcode_or_linker(self, reads, transposon_seq):
        # Instance without sequences, should extract genomic.
        extractor = BasicGenomicExtractor(
            transposon_seq, barcode_sequences=None,
            barcode_map=None, linker_sequence=None)

        genomic, barcode = extractor.extract_read(reads[3])
        assert genomic.sequence == 'GCGTCTG'
        assert barcode is None

    def test_exact_no_barcode_or_linker_neg(self, reads, transposon_seq,
                                            barcode_seqs, linker_seq):
        # Instance with barcodes and linker sequence, should return None.
        extractor = BasicGenomicExtractor(
            transposon_seq, barcode_sequences=barcode_seqs,
            barcode_map=None, linker_sequence=linker_seq)

        res = extractor.extract_read(reads[3])
        assert res is None


# -------- Identifier tests -------

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


# noinspection PyShadowingNames
# noinspection PyMethodMayBeStatic
class TestInsertionIdentifier(object):

    def test_group_alignments_by_position(self, test_data):
        alignments, _ = test_data

        identifier = InsertionIdentifier()
        groups = list(identifier._group_alignments_by_position(alignments))

        # Should identify a total of two groups.
        assert len(groups) == 2

        # Check group 1 for position and membership.
        (pos, strand, bc), alignments = groups[0]
        assert pos == 20
        assert strand == 1
        assert np.isnan(bc)
        assert len(alignments) == 4

    def test_group_alignments_by_position_bc(self, test_data):
        alignments, bc_map = test_data

        identifier = InsertionIdentifier()
        groups = list(identifier._group_alignments_by_position(
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

    def test_group_alignments_by_position_stranded(self, test_data_stranded):
        alignments = test_data_stranded

        identifier = InsertionIdentifier()
        groups = list(identifier._group_alignments_by_position(alignments))

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

    def test_group_alignments_by_position_unordered(self, test_data_unordered):
        with pytest.raises(ValueError):
            identifier = InsertionIdentifier()
            list(identifier._group_alignments_by_position(test_data_unordered))
