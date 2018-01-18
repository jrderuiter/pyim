from genopandas import GenomicDataFrame
import pytest

from pyim.annotate.util import filter_matches, annotate_matches
from pyim.model import InsertionSet


@pytest.fixture
def insertions():
    """Example set of insertions."""
    file_path = pytest.helpers.data_path('tagmap/subset.txt')
    return InsertionSet.from_csv(file_path, sep='\t')


@pytest.fixture
def matches():
    """Example matches for above insertions."""
    file_path = pytest.helpers.data_path('tagmap/subset.ann.txt')
    insertions = InsertionSet.from_csv(file_path, sep='\t')
    return insertions[['id', 'gene_id']].drop_duplicates()


@pytest.fixture
def matches_ann():
    """Example matches for insertions, with annotations."""
    file_path = pytest.helpers.data_path('tagmap/subset.ann.txt')
    insertions = InsertionSet.from_csv(file_path, sep='\t')

    cols = ['id', 'gene_id', 'gene_name', 'gene_distance', 'gene_orientation']
    return insertions[cols].drop_duplicates()


@pytest.fixture(scope='module')
def genes():
    """Example reference gene set."""
    file_path = pytest.helpers.data_path('tagmap/reference.gtf.gz')
    return GenomicDataFrame.from_gtf(
        file_path, filter_=lambda rec: rec['feature'] == 'gene')


class TestFilterMatches(object):
    """Tests for filter_matches function."""

    def test_closest(self, matches_ann):
        """Tests filtering for closest gene."""
        filtered = filter_matches(matches_ann, closest=True)

        # Check gene annotations.
        annotations = {
            id_: set(grp['gene_name'])
            for id_, grp in filtered.groupby('id')
        }

        assert annotations == {
            'test.INS_2': {'Capn9'},
            'test.INS_5': {'Trps1'},
            'test.INS_6': {'Slc16a9'},
            'test.INS_7': {'Ppp1r12a'},
            'test.INS_4': {'Arid1a'},
            'test.INS_3': {'Fgfr2'},
            'test.INS_1': {'Tanc1'}
        }

    def test_blacklist(self, matches_ann):
        """Tests filtering with blacklist."""

        filtered = filter_matches(
            matches_ann,
            closest=False,
            blacklist={'Slc16a9', 'Capn9', 'Trps1'},
            blacklist_col='gene_name')

        # Check gene annotations.
        annotations = {
            id_: set(grp['gene_name'])
            for id_, grp in
            filtered.dropna(subset=['gene_name'])
                    .groupby('id')
        }  # yapf: disable

        assert annotations == {
            'test.INS_2': {'Agt'},
            'test.INS_6': {'Gm7097', 'Gm9877'},
            'test.INS_7': {'Ppp1r12a'},
            'test.INS_4': {'Arid1a'},
            'test.INS_3': {'Fgfr2'},
            'test.INS_1': {'Tanc1'}
        }

    def test_blacklist_closest(self, matches_ann):
        """Tests filtering for closest genes with blacklist."""

        filtered = filter_matches(
            matches_ann,
            closest=True,
            blacklist={'Slc16a9', 'Capn9', 'Trps1'},
            blacklist_col='gene_name')

        # Check gene annotations.
        annotations = {
            id_: set(grp['gene_name'])
            for id_, grp in
            filtered.dropna(subset=['gene_name'])
                    .groupby('id')
        }  # yapf: disable

        assert annotations == {
            'test.INS_7': {'Ppp1r12a'},
            'test.INS_4': {'Arid1a'},
            'test.INS_3': {'Fgfr2'},
            'test.INS_1': {'Tanc1'}
        }


class TestAnnotateMatches(object):
    """Tests for annotate_matches function."""

    def test_example(self, matches, insertions, genes):
        """Tests annotation of matches."""

        annotated = annotate_matches(
            matches, insertions=insertions, genes=genes)

        assert 'gene_name' in annotated
        assert 'gene_id' in annotated
        assert 'gene_name' in annotated

        assert len(annotated) == len(matches)
