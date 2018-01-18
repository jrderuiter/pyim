from genopandas import GenomicDataFrame
import pytest

from pyim.annotate.annotators import rbm
from pyim.model import InsertionSet


class TestRbmAnnotator(object):
    """Tests for the RbmAnnotator class."""

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
        annotator = rbm.RbmAnnotator(preset='SB')
        matches = annotator.match(insertions, genes)

        # Cehck annotations.
        annotations = {
            id_: set(zip(grp['gene_id'], grp['window']))
            for id_, grp in matches.groupby('id')
        }

        assert annotations == {
            'test.INS_3': {('ENSMUSG00000030849', 'is')},
            'test.INS_6': {('ENSMUSG00000037762', 'ia')},
            'test.INS_2': {('ENSMUSG00000031980', 'ds'),
                           ('ENSMUSG00000031981', 'ia')},
            'test.INS_1': {('ENSMUSG00000035168', 'ia')},
            'test.INS_4': {('ENSMUSG00000007880', 'ia')},
            'test.INS_7': {('ENSMUSG00000019907', 'ia')},
            'test.INS_5': {('ENSMUSG00000038679', 'ia')}
        }  # yapf: disable

    def test_annotate(self, insertions, genes):
        """Tests annotate method."""

        # Annotate insertions.
        annotator = rbm.RbmAnnotator(preset='SB')
        annotated = annotator.annotate(insertions, genes)

        # Check gene annotations.
        annotations = {
            id_: set(zip(grp['gene_name'], grp['window']))
            for id_, grp in annotated.groupby('id')
        }

        assert annotations == {
            'test.INS_3': {('Fgfr2', 'is')},
            'test.INS_7': {('Ppp1r12a', 'ia')},
            'test.INS_2': {('Agt', 'ds'), ('Capn9', 'ia')},
            'test.INS_1': {('Tanc1', 'ia')},
            'test.INS_6': {('Slc16a9', 'ia')},
            'test.INS_4': {('Arid1a', 'ia')},
            'test.INS_5': {('Trps1', 'ia')}
        }
