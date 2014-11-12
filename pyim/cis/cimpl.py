__author__ = 'Julian'

import pandas.rpy.common as com

from rpy2 import robjects
from rpy2.robjects.packages import importr


def annotate_with_cis(ins_frame, genome='mm10', chromosomes=None,
                      scales=30000, n_iter=10, alpha=0.05, verbose=False):
    cimpl = importr('cimpl')

    ## Convert frame to cimpl format in R.
    ins_cimpl = _convert_to_cimpl_frame(ins_frame)
    ins_cimpl_r = com.convert_to_r_dataframe(ins_cimpl)

    ## Run the CIMPL analysis.
    if chromosomes is None:
        chromosomes = list(ins_cimpl['chr'].unique())
    chromosomes = robjects.vectors.StrVector(chromosomes)

    if type(scales) == list:
        scales = robjects.vectors.IntVector(scales)

    genome_obj = _load_genome(genome)

    cimpl_result = cimpl.doCimplAnalysis(
        ins_cimpl_r, scales=scales, n_iterations=n_iter, system='SB',
        BSgenome=genome_obj, chromosomes=chromosomes, verbose=verbose)

    ## Determine which insertions belong to CISs.
    cis_ins_ids = _get_cis_insertions(cimpl_result, alpha=alpha, mul_test=True)

    ## Create annotated copy of our insertion data frame.
    ins_copy = ins_frame.copy(deep=True)
    ins_copy['cis'] = ins_copy['id'].isin(set(cis_ins_ids))

    return ins_copy


def _convert_to_cimpl_frame(ins_frame):
    ## Extract and rename required columns.
    cimpl = ins_frame.ix[:,['id', 'seqname', 'location', 'sample']]
    cimpl.columns = ['id', 'chr', 'location', 'sampleID']

    ## Add 'chr' prefix to the chromosome names.
    cimpl['chr'] = ['chr' + c for c in cimpl['chr']]

    return cimpl


def _load_genome(genome):
    if genome == 'mm10':
        bs_genome = importr('BSgenome.Mmusculus.UCSC.mm10')
        genome_obj = bs_genome.Mmusculus
    else:
        raise ValueError('Unknown genome reference {}'.format(genome))

    return genome_obj


def _get_cis_insertions(cimpl_result, alpha=0.05, mul_test=True):
    cimpl = importr('cimpl')

    ## Identify CISs using CIMPL.
    cis = cimpl.getCISs(cimpl_result, alpha=alpha, mul_test=mul_test)

    ## Extract scale information from CIMPL result.
    cis_matrix_r = cimpl.getCISMatrix(cimpl_result, cis)
    cis_matrix = com.convert_robj(cis_matrix_r)
    cis_matrix_scales = cis_matrix.ix[:,[c.startswith('X') for c in cis_matrix.columns]]

    ## Select insertions that are attributed to at least
    ## one CIS by CIMPL on any scale.
    cis_mask = (cis_matrix_scales != '').apply(any, axis=1)
    cis_ins_ids = cis_matrix.ix[cis_mask, 'id']

    return list(cis_ins_ids)


##ins = pandas.read_csv('/Users/Julian/Desktop/Spontaenous.insertions.ulp2.txt', sep='\t')
##ins_annot = annotate_with_cis(ins, chromosomes=['chr10'])
