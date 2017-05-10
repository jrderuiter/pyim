import pytest

from pyim.annotate.annotators.rbm import RbmAnnotator, RbmCisAnnotator

# pylint: disable=redefined-outer-name


class TestWindowAnnotator(object):
    def test_basic(self, insertions, gtf_path):
        annotator = RbmAnnotator(gtf_path, preset='SB', verbose=False)
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


class TestWindowCisAnnotator(object):
    def test_basic(self, cis_insertions, cis_sites, gtf_path):
        annotator = RbmCisAnnotator(
            gtf_path, preset='SB', cis_sites=cis_sites, verbose=False)

        annotated = {ins.id: ins for ins in annotator.annotate(cis_insertions)}

        assert annotated['INS1'].metadata['cis_id'] == 'CIS1'
        assert annotated['INS1'].metadata['gene_name'] == 'Trp53bp2'

        assert annotated['INS2'].metadata['cis_id'] == 'CIS2'
        assert 'gene_name' not in annotated['INS2'].metadata
