from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
from builtins import (ascii, bytes, chr, dict, filter, hex, input,
                      int, map, next, oct, open, pow, range, round,
                      str, super, zip)

from collections import namedtuple

import pytest

import pysam
from skbio import DNASequence

from pyim.pipelines.shear_splink import \
    ShearSplinkExtractor, ShearSplinkStatus, ShearSplinkIdentifier


# --- Extractor --- #

@pytest.fixture(scope='module')
def reads():
    return [
        DNASequence('AAAT' + 'ACTG' + 'GCGTCTG' + 'CCCC'),  # BC01, tr, linker
        DNASequence('AAAA' + 'ACTG' + 'GCGTCTG' + 'CCCC'),  # BC02, tr, linker
        DNASequence('AAAT' + 'ACTG' + 'GCGTCTG'),           # No linker
        DNASequence('ACTG' + 'GCGTCTG' + 'CCCC'),           # No barcode
        DNASequence('AAAA' + 'GCGTCTG' + 'CCCC')            # No transposon
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
class TestShearSplinkExtractor(object):

    def test_forward(self, reads, transposon_seq, barcode_seqs, linker_seq):
        extractor = ShearSplinkExtractor(
            transposon_sequence=transposon_seq,
            barcode_sequences=barcode_seqs,
            linker_sequence=linker_seq)

        (genomic, barcode), status = extractor.extract_read(reads[0])
        assert genomic.sequence == 'GCGTCTG'
        assert barcode == 'BC01'
        assert status == ShearSplinkStatus.proper_read

        (genomic, barcode), status = extractor.extract_read(reads[1])
        assert genomic.sequence == 'GCGTCTG'
        assert barcode == 'BC02'
        assert status == ShearSplinkStatus.proper_read

    def test_reverse(self, reads, transposon_seq, barcode_seqs, linker_seq):
        extractor = ShearSplinkExtractor(
            transposon_sequence=transposon_seq,
            barcode_sequences=barcode_seqs,
            linker_sequence=linker_seq)

        res, status = extractor.extract_read(reads[0].reverse_complement())
        assert res is not None
        assert status == ShearSplinkStatus.proper_read

        genomic, barcode = res
        assert genomic.sequence == 'GCGTCTG'
        assert barcode == 'BC01'

    def test_missing_linker(self, reads, transposon_seq,
                            barcode_seqs, linker_seq):
        extractor = ShearSplinkExtractor(
            transposon_sequence=transposon_seq,
            barcode_sequences=barcode_seqs,
            linker_sequence=linker_seq)

        res, status = extractor.extract_read(reads[2])
        assert res is None
        assert status == ShearSplinkStatus.no_linker

    def test_missing_barcode(self, reads, transposon_seq,
                             barcode_seqs, linker_seq):
        extractor = ShearSplinkExtractor(
            transposon_sequence=transposon_seq,
            barcode_sequences=barcode_seqs,
            linker_sequence=linker_seq)

        res, status = extractor.extract_read(reads[3])
        assert res is None
        assert status == ShearSplinkStatus.no_barcode

    def test_missing_transposon(self, reads, transposon_seq,
                                barcode_seqs, linker_seq):
        extractor = ShearSplinkExtractor(
            transposon_sequence=transposon_seq,
            barcode_sequences=barcode_seqs,
            linker_sequence=linker_seq)

        res, status = extractor.extract_read(reads[4])
        assert res is None
        assert status == ShearSplinkStatus.no_transposon

# --- Identifier --- #


AlignedSegment = namedtuple(
    'AlignedSegment',
    ['query_name', 'reference_id', 'reference_start',
     'reference_end', 'is_reverse', 'mapping_quality'])


class AlignmentFile(object):

    def __init__(self, *args, **kwargs):
        pass

    def fetch(self, reference=None):
        raise NotImplementedError()

    @property
    def references(self):
        raise NotImplementedError()


@pytest.fixture()
def patch_alignment_file(monkeypatch):
    alignments = {'1': [
        AlignedSegment('READ0', '1', 1, 21, True, 40),    # 1st group on -1.
        AlignedSegment('READ1', '1', 20, 70, False, 40),  # Start of 2nd group.
        AlignedSegment('READ2', '1', 20, 70, False, 40),
        AlignedSegment('READ3', '1', 20, 69, False, 40),
        AlignedSegment('READ5', '1', 21, 68, False, 40),  # Within merge dist.
        AlignedSegment('READ6', '1', 32, 68, False, 30)   # Outside merge dist.
    ]}

    # Monkeypatch mock class for fixture.
    monkeypatch.setattr(AlignmentFile, 'fetch',
                        lambda self, reference: alignments[reference])
    monkeypatch.setattr(AlignmentFile, 'references', list(alignments.keys()))

    # Monkeypatch class into pysam.
    monkeypatch.setattr(pysam, 'AlignmentFile', AlignmentFile)


# noinspection PyShadowingNames
# noinspection PyMethodMayBeStatic
class TestLamPcrIdentifier(object):

    def test_identify(self, patch_alignment_file):
        identifier = ShearSplinkIdentifier(merge_distance=10, min_mapq=0)
        insertions = identifier.identify('dummy.bam')

        assert len(insertions) == 3

    def test_identify_mapq(self, patch_alignment_file):
        identifier = ShearSplinkIdentifier(merge_distance=10, min_mapq=37)
        insertions = identifier.identify('dummy.bam')

        # Should find two insertions, due to lower mapq of READ6.
        assert len(insertions) == 2

    def test_identify_large_merge(self, patch_alignment_file):
        identifier = ShearSplinkIdentifier(merge_distance=20, min_mapq=0)
        insertions = identifier.identify('dummy.bam')

        # Should find two insertions, due to larger merge_dist.
        assert len(insertions) == 2

    def test_identify_no_merge(self, patch_alignment_file):
        identifier = ShearSplinkIdentifier(merge_distance=0, min_mapq=0)
        insertions = identifier.identify('dummy.bam')

        # Should find four insertions, due to no merge_dist.
        assert len(insertions) == 4
