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
        identifier = LamPcrIdentifier(merge_distance=10, min_mapq=0)
        insertions = identifier.identify('dummy.bam')

        assert len(insertions) == 3

    def test_identify_mapq(self, patch_alignment_file):
        identifier = LamPcrIdentifier(merge_distance=10, min_mapq=37)
        insertions = identifier.identify('dummy.bam')

        # Should find two insertions, due to lower mapq of READ6.
        assert len(insertions) == 2

    def test_identify_large_merge(self, patch_alignment_file):
        identifier = LamPcrIdentifier(merge_distance=20, min_mapq=0)
        insertions = identifier.identify('dummy.bam')

        # Should find two insertions, due to larger merge_dist.
        assert len(insertions) == 2

    def test_identify_no_merge(self, patch_alignment_file):
        identifier = LamPcrIdentifier(merge_distance=0, min_mapq=0)
        insertions = identifier.identify('dummy.bam')

        # Should find four insertions, due to no merge_dist.
        assert len(insertions) == 4
