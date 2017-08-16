import pytest

from pyim.annotate.annotators.rbm import RbmAnnotator

# pylint: disable=redefined-outer-name


class TestRbmAnnotator(object):
    """Unit tests for the RbmAnnotator."""

    def test_basic(self, insertions, genes):
        """Tests a simple example."""

        annotator = RbmAnnotator(genes, preset='SB')
        annotated = list(annotator.annotate(insertions))

        metadata1 = annotated[0].metadata
        assert metadata1['gene_name'] == 'Trp53bp2'
        assert metadata1['gene_id'] == 'ENSMUSG00000026510'
        assert metadata1['gene_distance'] == -1000
        assert metadata1['gene_orientation'] == 'sense'

        metadata2 = annotated[1].metadata
        assert metadata2['gene_name'] == 'Myh9'
        assert metadata2['gene_id'] == 'ENSMUSG00000022443'
        assert metadata2['gene_distance'] == 2000
        assert metadata2['gene_orientation'] == 'antisense'

        assert 'gene_name' not in annotated[2].metadata

    def test_select_closest(self, insertions, genes):
        """Test with selection for closest gene."""

        annotator = RbmAnnotator(genes, preset='SB', closest=True)
        annotated = list(annotator.annotate(insertions))

        assert len(annotated) == 3
        assert annotated[0].metadata['gene_name'] == 'Trp53bp2'
