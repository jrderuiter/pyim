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

GENOMIC_SEQ1 = DNASequence('CACTGGCCACGCGAAGGTGC')
GENOMIC_SEQ2 = DNASequence('GACCACTGGCCACGCGAAGG').reverse_complement()
GENOMIC_SEQ3 = DNASequence('CGTTGGTCACTCTACCCACA')


@pytest.fixture(scope='module')
def transposon_seq():
    return DNASequence('TTTG', id='transposon')


@pytest.fixture(scope='module')
def barcode_seqs():
    return [DNASequence('AAAT', 'BC01'),
            DNASequence('AAAA', 'BC02')]


@pytest.fixture(scope='module')
def linker_seq():
    return DNASequence('CCCG', id='linker')


# noinspection PyShadowingNames
# noinspection PyMethodMayBeStatic
class TestShearSplinkExtractor(object):

    def test_forward(self, transposon_seq, barcode_seqs, linker_seq):
        read = DNASequence(str(barcode_seqs[0]) + str(transposon_seq) +
                           str(GENOMIC_SEQ1) + str(linker_seq))

        extractor = ShearSplinkExtractor(
            transposon_sequence=transposon_seq,
            barcode_sequences=barcode_seqs,
            linker_sequence=linker_seq)

        result, status = extractor.extract_read(read)
        assert status == ShearSplinkStatus.proper_read
        assert result is not None

        genomic, barcode = result
        assert genomic.sequence == GENOMIC_SEQ1.sequence
        assert barcode == barcode_seqs[0].id

    def test_reverse(self, transposon_seq,
                     barcode_seqs, linker_seq):
        read = DNASequence(str(barcode_seqs[0]) + str(transposon_seq) +
                           str(GENOMIC_SEQ1) + str(linker_seq))
        read = read.reverse_complement()

        extractor = ShearSplinkExtractor(
            transposon_sequence=transposon_seq,
            barcode_sequences=barcode_seqs,
            linker_sequence=linker_seq)

        res, status = extractor.extract_read(read)
        assert res is not None
        assert status == ShearSplinkStatus.proper_read

        genomic, barcode = res
        assert genomic.sequence == GENOMIC_SEQ1.sequence
        assert barcode == barcode_seqs[0].id

    def test_missing_linker(self, transposon_seq, barcode_seqs, linker_seq):
        read = DNASequence(str(barcode_seqs[0]) + str(transposon_seq) +
                           str(GENOMIC_SEQ1))

        extractor = ShearSplinkExtractor(
            transposon_sequence=transposon_seq,
            barcode_sequences=barcode_seqs,
            linker_sequence=linker_seq)

        res, status = extractor.extract_read(read)
        assert res is None
        assert status == ShearSplinkStatus.no_linker

    def test_missing_barcode(self, transposon_seq,
                             barcode_seqs, linker_seq):
        read = DNASequence(str(transposon_seq) + str(GENOMIC_SEQ1) +
                           str(linker_seq))

        extractor = ShearSplinkExtractor(
            transposon_sequence=transposon_seq,
            barcode_sequences=barcode_seqs,
            linker_sequence=linker_seq)

        res, status = extractor.extract_read(read)
        assert res is None
        assert status == ShearSplinkStatus.no_barcode

    def test_missing_transposon(self, transposon_seq,
                                barcode_seqs, linker_seq):
        read = DNASequence(str(barcode_seqs[0]) + str(GENOMIC_SEQ1) +
                           str(linker_seq))

        extractor = ShearSplinkExtractor(
            transposon_sequence=transposon_seq,
            barcode_sequences=barcode_seqs,
            linker_sequence=linker_seq)

        res, status = extractor.extract_read(read)
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


@pytest.fixture()
def aln_barcode_map():
    return {
        'READ0': 'Sample1',
        'READ1': 'Sample1',
        'READ2': 'Sample1',
        'READ3': 'Sample2',
        'READ4': 'Sample2',
        'READ5': 'Sample2',
        'READ6': 'Sample2',
    }


# noinspection PyShadowingNames
# noinspection PyMethodMayBeStatic
class TestShearSplinkIdentifier(object):

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

    def test_identify_barcodes(self, patch_alignment_file, aln_barcode_map):
        identifier = ShearSplinkIdentifier(merge_distance=10, min_mapq=0)
        insertions = identifier.identify(
            'dummy.bam', barcode_map=aln_barcode_map)

        # Should find four insertions, as reads (1,2) and 3
        # now belong to different samples.
        assert len(insertions) == 4

        # Test insertion membership.
        assert insertions['sample'].iloc[0] == 'Sample1'
        assert insertions['sample'].iloc[1] == 'Sample2'
        assert insertions['sample'].iloc[2] == 'Sample1'
        assert insertions['sample'].iloc[3] == 'Sample2'
