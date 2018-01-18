from genopandas import GenomicDataFrame
import pytest

from pyim.annotate.annotators import window
from pyim.model import InsertionSet


class TestWindowAnnotator(object):
    """Tests for the WindowAnnotator class."""

    @pytest.fixture
    def insertions(self):
        """Example set of insertions."""
        file_path = pytest.helpers.data_path('tagmap/subset.txt')
        return InsertionSet.from_csv(file_path, sep='\t')

    @pytest.fixture(scope='module')
    def genes(self):
        """Example reference gene set."""
        file_path = pytest.helpers.data_path('tagmap/reference.gtf.gz')
        return GenomicDataFrame.from_gtf(
            file_path, filter_=lambda rec: rec['feature'] == 'gene')

    def test_match(self, insertions, genes):
        """Tests match method."""

        # Annotate insertions.
        annotator = window.WindowAnnotator.from_window_size(20000)
        matches = annotator.match(insertions, genes)

        assert list(matches.columns) == ['id', 'gene_id', 'window']

        # Check gene annotations.
        annotations = {
            id_: set(grp['gene_id'])
            for id_, grp in matches.groupby('id')
        }

        assert annotations == {
            'test.INS_7': {'ENSMUSG00000019907'},
            'test.INS_3': {'ENSMUSG00000030849'},
            'test.INS_2': {'ENSMUSG00000031980', 'ENSMUSG00000031981'},
            'test.INS_6':
            {'ENSMUSG00000052426', 'ENSMUSG00000037762', 'ENSMUSG00000089835'},
            'test.INS_1': {'ENSMUSG00000035168'},
            'test.INS_4': {'ENSMUSG00000007880'},
            'test.INS_5': {'ENSMUSG00000038679'}
        }

    def test_annotate(self, insertions, genes):
        """Tests annotate method."""

        # Annotate insertions.
        annotator = window.WindowAnnotator.from_window_size(20000)
        annotated = annotator.annotate(insertions, genes=genes)

        # Check returned type.
        assert isinstance(annotated, InsertionSet)

        # Check gene annotations.
        annotations = {
            id_: set(grp['gene_id'])
            for id_, grp in annotated.groupby('id')
        }

        assert annotations == {
            'test.INS_7': {'ENSMUSG00000019907'},
            'test.INS_3': {'ENSMUSG00000030849'},
            'test.INS_2': {'ENSMUSG00000031980', 'ENSMUSG00000031981'},
            'test.INS_6':
            {'ENSMUSG00000052426', 'ENSMUSG00000037762', 'ENSMUSG00000089835'},
            'test.INS_1': {'ENSMUSG00000035168'},
            'test.INS_4': {'ENSMUSG00000007880'},
            'test.INS_5': {'ENSMUSG00000038679'}
        }

        # Check a few values.
        ann_ind = annotated.set_index(['id', 'gene_name'])

        assert ann_ind.loc[('test.INS_2', 'Agt'), 'gene_distance'] == 13709
        assert ann_ind.loc[('test.INS_2', 'Agt'), 'gene_orientation'] == 1

        assert ann_ind.loc[('test.INS_2', 'Capn9'), 'gene_distance'] == 0
        assert ann_ind.loc[('test.INS_2', 'Capn9'), 'gene_orientation'] == -1

    def test_annotate_closest(self, insertions, genes):
        """Tests annotate with closest == True."""

        # Annotate insertions.
        annotator = window.WindowAnnotator.from_window_size(20000)
        annotated = annotator.annotate(insertions, genes=genes, closest=True)

        # Check gene annotations.
        annotations = {
            id_: set(grp['gene_name'])
            for id_, grp in annotated.groupby('id')
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

    def test_annotate_blacklist(self, insertions, genes):
        """Tests annotate with blacklist."""

        # Annotate insertions.
        annotator = window.WindowAnnotator.from_window_size(20000)

        annotated = annotator.annotate(
            insertions,
            genes=genes,
            closest=False,
            blacklist={'Slc16a9', 'Capn9', 'Trps1'},
            blacklist_col='gene_name')

        # Check gene annotations.
        annotations = {
            id_: set(grp['gene_name'])
            for id_, grp in
            annotated.dropna(subset=['gene_name'])
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

    def test_annotate_closest_blacklist(self, insertions, genes):
        """Tests annotate with closest == True and blacklist.

        Should drop INS_2 and INS_6 entirely, as blacklisting is
        done after selecting for the closest gene.
        """

        # Annotate insertions.
        annotator = window.WindowAnnotator.from_window_size(20000)

        annotated = annotator.annotate(
            insertions,
            genes=genes,
            closest=True,
            blacklist={'Slc16a9', 'Capn9', 'Trps1'},
            blacklist_col='gene_name')

        # Check gene annotations.
        annotations = {
            id_: set(grp['gene_name'])
            for id_, grp in
            annotated.dropna(subset=['gene_name'])
                    .groupby('id')
        }  # yapf: disable

        assert annotations == {
            'test.INS_7': {'Ppp1r12a'},
            'test.INS_4': {'Arid1a'},
            'test.INS_3': {'Fgfr2'},
            'test.INS_1': {'Tanc1'}
        }
