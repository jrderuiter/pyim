from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
from builtins import (ascii, bytes, chr, dict, filter, hex, input,
                      int, map, next, oct, open, pow, range, round,
                      str, super, zip)

import pytest
from skbio import DNASequence

from pyim.pipelines.shear_splink import ShearSplinkExtractor, ShearSplinkStatus


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
