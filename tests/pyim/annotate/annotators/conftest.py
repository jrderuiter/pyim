import pytest
from pyim.util.frozendict import frozendict

from pyim.model import Insertion, CisSite


@pytest.fixture
def gtf_path():
    return pytest.helpers.data_path('reference.gtf')


@pytest.fixture
def insertions():
    # Trp53bp2 location: 1: 182,409,172-182,462,432.
    # Myh9 location: 15: 77,760,587-77,842,175.

    return [
        # 1000 bp upstream of Trp53bp2.
        Insertion(id='INS1', chromosome='1', position=182408172,
                  strand=1, support=2, metadata=frozendict()),
        # 2000 bp downstream of Myh9.
        Insertion(id='INS2', chromosome='15', position=77758587,
                  strand=1, support=2, metadata=frozendict()),
        # Different chromosome.
        Insertion(id='INS3', chromosome='4', position=77843175,
                  strand=1, support=2, metadata=frozendict())
    ] # yapf: disable


@pytest.fixture
def cis_insertions():
    return [
        # 1000 bp upstream of Trp53bp2.
        Insertion(id='INS1', chromosome='1', position=182408172, strand=1,
                support=2, metadata=frozendict({'cis_id': 'CIS1'})),
        # Different chromosome.
        Insertion(id='INS2', chromosome='4', position=77843175, strand=1,
                  support=2, metadata=frozendict({'cis_id': 'CIS2'}))
    ] # yapf: disable


@pytest.fixture
def cis_sites():
    return [
        CisSite(id='CIS1', chromosome='1', position=182408172,
                strand=1, metadata=frozendict()),
        CisSite(id='CIS2', chromosome='4', position=132408091,
                strand=1, metadata=frozendict())
    ] # yapf: disable
