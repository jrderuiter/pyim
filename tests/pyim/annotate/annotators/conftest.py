from pathlib import Path

import pytest

from pyim.model import Insertion, CisSite
from pyim.util.frozendict import frozendict
from pyim.util.pandas import GenomicDataFrame


@pytest.fixture(scope='session')
def gtf_path():
    """Path to example GTF file."""
    return Path(str(pytest.helpers.data_path('reference.gtf.gz')))


@pytest.fixture(scope='session')
def genes(gtf_path):
    """Genes from example GTF file."""
    genes = GenomicDataFrame.from_gtf(gtf_path, feature='gene')
    genes['strand'] = genes['strand'].map({'+': 1, '-': -1})
    return genes


@pytest.fixture(scope='session')
def insertions():
    """Example insertions."""

    # Trp53bp2 location: 1: 182,409,172-182,462,432.
    # Myh9 location: 15: 77,760,587-77,842,175.

    return [
        # 1000 bp upstream of Trp53bp2.
        Insertion(id='INS1', chromosome='1', position=182408171,
                  strand=1, support=2, sample='s1', metadata=frozendict()),
        # 2000 bp downstream of Myh9.
        Insertion(id='INS2', chromosome='15', position=77758586,
                  strand=1, support=2, sample='s1', metadata=frozendict()),
        # Different chromosome.
        Insertion(id='INS3', chromosome='4', position=77843175,
                  strand=1, support=2, sample='s1', metadata=frozendict())
    ] # yapf: disable


@pytest.fixture(scope='session')
def cis_insertions():
    """Example insertions with CIS annotations."""

    return [
        # 1000 bp upstream of Trp53bp2.
        Insertion(id='INS1', chromosome='1', position=182408172,
                  strand=1, support=2, sample='s1',
                  metadata=frozendict({'cis_id': 'CIS1'})),
        # Different chromosome.
        Insertion(id='INS2', chromosome='4', position=77843175,
                  strand=1, support=2, sample='s1',
                  metadata=frozendict({'cis_id': 'CIS2'}))
    ] # yapf: disable


@pytest.fixture(scope='session')
def cis_sites():
    """Example CIS sites."""

    return [
        CisSite(id='CIS1', chromosome='1', position=182408172,
                strand=1, metadata=frozendict()),
        CisSite(id='CIS2', chromosome='4', position=132408091,
                strand=1, metadata=frozendict())
    ] # yapf: disable
