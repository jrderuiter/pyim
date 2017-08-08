from pyim.annotate.annotators import CisAnnotator, WindowAnnotator


class TestCisAnnotator(object):
    def test_basic(self, cis_insertions, cis_sites, genes):

        annotator = CisAnnotator(
            annotator=WindowAnnotator.from_window_size(
                genes=genes, window_size=1000),
            genes=genes,
            cis_sites=cis_sites)

        annotated = {ins.id: ins for ins in annotator.annotate(cis_insertions)}

        assert annotated['INS1'].metadata['cis_id'] == 'CIS1'
        assert annotated['INS1'].metadata['gene_name'] == 'Trp53bp2'

        assert annotated['INS2'].metadata['cis_id'] == 'CIS2'
        assert 'gene_name' not in annotated['INS2'].metadata
