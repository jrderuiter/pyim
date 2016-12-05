from pyim.annotate.annotators.window import (Window, WindowAnnotator,
                                             WindowCisAnnotator)

# pylint: disable=redefined-outer-name


def build_window(size):
    return Window(
        -(size // 2), (size // 2),
        strand=None,
        name=None,
        strict_left=False,
        strict_right=False)


class TestWindowAnnotator(object):
    def test_basic(self, insertions, gtf_path):
        annotator = WindowAnnotator(
            gtf_path, [build_window(20000)], verbose=False)
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

    def test_small_window(self, insertions, gtf_path):
        annotator = WindowAnnotator(
            gtf_path, [build_window(2000)], verbose=False)
        annotated = list(annotator.annotate(insertions))

        assert annotated[0].metadata['gene_name'] == 'Trp53bp2'
        assert 'gene_name' not in annotated[1].metadata
        assert 'gene_name' not in annotated[2].metadata

    def test_blacklist(self, insertions, gtf_path):
        annotator = WindowAnnotator(
            gtf_path, [build_window(20000)],
            verbose=False,
            blacklist={'Trp53bp2'})
        annotated = list(annotator.annotate(insertions))

        assert 'gene_name' not in annotated[0].metadata

    def test_multiple(self, insertions, gtf_path):
        annotator = WindowAnnotator(
            gtf_path, [build_window(200000000)], verbose=False)
        annotated = list(annotator.annotate(insertions))

        assert len(annotated) == 4

        assert annotated[0].id == 'INS1'
        assert annotated[1].id == 'INS1'

        genes = {annotated[0].metadata['gene_name'],
                 annotated[1].metadata['gene_name']}
        assert genes == {'Ppp1r12b', 'Trp53bp2'}

    def test_select_closest(self, insertions, gtf_path):
        annotator = WindowAnnotator(
            gtf_path, [build_window(200000000)], closest=True, verbose=False)
        annotated = list(annotator.annotate(insertions))

        assert len(annotated) == 3
        assert annotated[0].metadata['gene_name'] == 'Trp53bp2'


class TestWindowCisAnnotator(object):
    def test_basic(self, cis_insertions, cis_sites, gtf_path):
        annotator = WindowCisAnnotator(
            gtf_path, [build_window(20000)],
            cis_sites=cis_sites,
            verbose=False)

        annotated = {ins.id: ins for ins in annotator.annotate(cis_insertions)}

        assert annotated['INS1'].metadata['cis_id'] == 'CIS1'
        assert annotated['INS1'].metadata['gene_name'] == 'Trp53bp2'

        assert annotated['INS2'].metadata['cis_id'] == 'CIS2'
        assert 'gene_name' not in annotated['INS2'].metadata
