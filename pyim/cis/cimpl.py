__author__ = 'Julian'


import readline  # Work-around for: ../lib/libreadline.so.6: undefined symbol: PC

from rpy2 import robjects
from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import importr

pandas2ri.activate()


def annotate_with_cis(ins_frame, scales, system=None, specificity_pattern=None, genome='mm10',
                      chromosomes=None, n_iter=10, alpha=0.05, verbose=True, threads=1):
    cimpl = importr('cimpl')

    # Convert frame to cimpl format in R.
    ins_cimpl = _convert_to_cimpl_frame(ins_frame)

    # Run the CIMPL analysis.
    if chromosomes is None:
        chromosomes = list(ins_cimpl['chr'].unique())
    chromosomes = robjects.vectors.StrVector(chromosomes)

    if type(scales) == list:
        scales = robjects.vectors.IntVector(scales)

    genome_obj = _load_genome(genome)

    # Actually run CIMPL, using either system or specificity pattern.
    if system is not None:
        cimpl_result = cimpl.doCimplAnalysis(
            ins_cimpl, scales=scales, n_iterations=n_iter, system=system,
            BSgenome=genome_obj, chromosomes=chromosomes, verbose=verbose, threads=threads)
    elif specificity_pattern is not None:
        cimpl_result = cimpl.doCimplAnalysis(
            ins_cimpl, scales=scales, n_iterations=n_iter,
            specificity_pattern=specificity_pattern, BSgenome=genome_obj,
            chromosomes=chromosomes, verbose=verbose, threads=threads)
    else:
        raise ValueError('Either the system or specificity pattern should be supplied.')

    # Determine which insertions belong to CISs.
    cis, cis_map = _get_cis(cimpl_result, alpha=alpha, mul_test=True)

    # Create annotated copy of our insertion data frame.
    cis_ins = ins_frame.ix[ins_frame['name'].isin(cis_map.index), :].copy()
    cis_ins['cis'] = [_concat_scales(cis_map.ix[n]) for n in cis_ins['name']]

    return cis_ins, cis


def _concat_scales(row):
    return ';'.join([sc for sc in row if type(sc) == str])


def _convert_to_cimpl_frame(ins_frame):
    # Extract and rename required columns.
    cimpl = ins_frame.ix[:, ['name', 'seqname', 'location', 'sample']]
    cimpl.columns = ['id', 'chr', 'location', 'sampleID']

    # Add 'chr' prefix to the chromosome names.
    cimpl['chr'] = ['chr' + c for c in cimpl['chr']]

    return cimpl


def _load_genome(genome):
    if genome == 'mm10':
        bs_genome = importr('BSgenome.Mmusculus.UCSC.mm10')
        genome_obj = bs_genome.Mmusculus
    else:
        raise ValueError('Unknown genome reference {}'.format(genome))

    return genome_obj


def _get_cis(cimpl_result, alpha=0.05, mul_test=True):
    cimpl = importr('cimpl')

    # Get cis entries from CIMPL.
    cis = cimpl.getCISs(cimpl_result, alpha=alpha, mul_test=mul_test)
    cis = robjects.r.cbind(cis, name=robjects.r.rownames(cis))
    cis = pandas2ri.ri2py_dataframe(cis)

    lookup = dict(zip(map(str, range(1, cis.shape[0] + 1)), cis['name']))

    # Extract matrix mapping insertions to CISs by index.
    cis_matrix = cimpl.getCISMatrix(cimpl_result, cis)
    cis_matrix = pandas2ri.ri2py_dataframe(cis_matrix)

    # Convert to DF mapping insertions --> CISs by name.
    cis_map = cis_matrix.set_index('id')
    cis_map = cis_map.ix[:, [c.startswith('X') for c in cis_map.columns]]

    cis_map = cis_map.apply(lambda x: x.map(lookup))
    cis_map.dropna(how='all', inplace=True)

    return cis, cis_map
