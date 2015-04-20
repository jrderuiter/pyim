from collections import namedtuple

import pytest

import pysam
from pyim.pipelines.lam_pcr import LamPcrIdentifier


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
    alignments = {0: [
        AlignedSegment('READ1', 0, 20, 70, False, 40),
        AlignedSegment('READ2', 0, 20, 70, False, 40),
        AlignedSegment('READ3', 0, 20, 69, False, 40),
        AlignedSegment('READ4', 0, 20, 68, False, 40),
        AlignedSegment('READ5', 0, 21, 68, False, 40)
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
        identifier = LamPcrIdentifier(merge_distance=10)
        #print(list(identifier.identify('alignment.bam')))