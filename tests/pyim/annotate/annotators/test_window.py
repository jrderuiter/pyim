from pyim.annotate.annotators.window import Window, WindowAnnotator

# pylint: disable=redefined-outer-name


class TestWindowAnnotator(object):
    """Unit tests for annotation using the WindowAnnotator."""

    def test_basic(self, insertions, genes):
        """Test simple example of annotating three insertions.

        The first two insertions should be annotated with Trp53bp2/Myh9
        respectively. The third should not be annotated with any gene,
        as this insertion occurs on a chromosome that does not contain
        any genes.

        Note that this should NOT raise an error for the non-existing
        chromosome.
        """

        # Annotate insertions.
        annotator = WindowAnnotator.from_window_size(
            genes=genes, window_size=20000)
        annotated = list(annotator.annotate(insertions))

        # Verify annotations.
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

        # Check that last insertion was not annotated.
        assert 'gene_name' not in annotated[2].metadata

    def test_small_window(self, insertions, genes):
        """Test example with smaller window. Should not annotate Myh9."""

        annotator = WindowAnnotator.from_window_size(genes, window_size=1500)
        annotated = list(annotator.annotate(insertions))

        assert annotated[0].metadata['gene_name'] == 'Trp53bp2'
        assert 'gene_name' not in annotated[1].metadata
        assert 'gene_name' not in annotated[2].metadata

    def test_blacklist(self, insertions, genes):
        """Test with blacklisted gene. Should not annotate Trp53bp2."""

        annotator = WindowAnnotator.from_window_size(
            genes, window_size=20000, blacklist={'ENSMUSG00000026510'})
        annotated = list(annotator.annotate(insertions))

        assert 'gene_name' not in annotated[0].metadata

    def test_multiple(self, insertions, genes):
        """Test with multiple gene annotations."""

        annotator = WindowAnnotator.from_window_size(
            genes, window_size=200000000)
        annotated = list(annotator.annotate(insertions))

        assert len(annotated) == 4

        assert annotated[0].id == 'INS1'
        assert annotated[1].id == 'INS1'

        genes = {
            annotated[0].metadata['gene_name'],
            annotated[1].metadata['gene_name']
        }
        assert genes == {'Ppp1r12b', 'Trp53bp2'}

    def test_select_closest(self, insertions, genes):
        """Test with selection for closest gene."""

        annotator = WindowAnnotator.from_window_size(
            genes, window_size=200000000, closest=True)
        annotated = list(annotator.annotate(insertions))

        assert len(annotated) == 3
        assert annotated[0].metadata['gene_name'] == 'Trp53bp2'
