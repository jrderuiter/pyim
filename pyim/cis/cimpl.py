__author__ = 'Julian'

import pandas as pd
from rpy2 import robjects

from tkgeno.util.rpy2 import importr, pandas_to_dataframe, dataframe_to_pandas


R_GENOMES = {
    'mm10': 'BSgenome.Mmusculus.UCSC.mm10'
}


def cimpl(insertions, scales, genome, system=None, pattern=None,
          lhc_method='none', iterations=1000, chromosomes=None,
          verbose=False, threads=1):
    # Fill in chromosomes from data if not specified.
    if chromosomes is None:
        chromosomes = list(insertions['seqname'].unique())

    # Determine if system or specific pattern was specified.
    if pattern is not None:
        extra_args = {'specificity_pattern': pattern}
    elif system is not None:
        extra_args = {'system': system}
    else:
        raise ValueError('Either system or specificity pattern '
                         'should be specified.')

    # Prepare chromosomes argument, adding 'chr' prefix and
    # converting to StrVector to pass to R.
    if not chromosomes[0].startswith('chr'):
        chromosomes = ['chr' + c for c in chromosomes]
    chromosomes = robjects.vectors.StrVector(chromosomes)

    # Convert scales to IntVector if supplied as list.
    if type(scales) == list:
        scales = robjects.vectors.IntVector(scales)

    # Load genome object from R.
    genome_obj = _load_genome(genome)

    # Run CIMPL!
    cimpl_r = importr('cimpl')
    cimpl_obj = cimpl_r.doCimplAnalysis(
        _convert_to_cimpl_dataframe(insertions),
        scales=scales, n_iterations=iterations,
        lhc_method=lhc_method, threads=threads, BSgenome=genome_obj,
        chromosomes=chromosomes, verbose=verbose, **extra_args)

    return cimpl_obj


def _convert_to_cimpl_dataframe(insertions):
    # Extract and rename required columns.
    cimpl_frame = insertions.ix[:, ['insertion_id', 'seqname',
                                    'location', 'sample']]
    cimpl_frame.columns = ['id', 'chr', 'location', 'sampleID']

    # Add 'chr' prefix to the chromosome names if needed.
    cimpl_frame['chr'] = _prefix_chromosomes(cimpl_frame['chr'])

    return pandas_to_dataframe(cimpl_frame)


def _prefix_chromosomes(series, prefix='chr'):
    # Add 'chr' prefix to the chromosome names if needed.
    if not series[0].startswith('chr'):
        series = series.map(lambda c: prefix + c)
    return series


def _load_genome(genome):
    # Lookup R package for genome.
    try:
        genome_pkg = R_GENOMES[genome]
    except KeyError:
        raise ValueError('Unsupported genome {}'.format(genome))

    # Import package and extract genome object.
    bs_genome = importr(genome_pkg)
    genome_obj = bs_genome.Mmusculus

    return genome_obj


def get_cis(cimpl_obj, alpha=0.05, mul_test=True):
    cimpl_r = importr('cimpl')
    cis_obj = cimpl_r.getCISs(cimpl_obj, alpha=alpha, mul_test=mul_test)

    # Convert cis to pandas and rename index.
    cis_frame = dataframe_to_pandas(cis_obj).reset_index()
    cis_frame.rename(columns={'index': 'cis_id',
                              'chromosome': 'seqname'}, inplace=True)

    # Convert columns to int types.
    for col in ['peak_location', 'start', 'end', 'width', 'n_insertions']:
        cis_frame[col] = cis_frame[col].astype(int)

    # Remove chr prefix from chromosomes.
    cis_frame['seqname'] = cis_frame['seqname'].str.replace('chr', '')

    # Reorder columns.
    cis_frame = cis_frame[['cis_id', 'seqname', 'start', 'end',
                           'scale', 'p_value', 'n_insertions',
                           'peak_location', 'peak_height', 'width']]

    return cis_frame


def get_cis_mapping(cimpl_obj, cis_frame):
    # Add cis_id as index to cis frame before passing to R,
    # ensures CIMPL uses cis id's instead of row indices.
    cis_frame = cis_frame.copy()
    cis_frame.set_index('cis_id', drop=False, inplace=True)
    cis_frame['chromosomes'] = _prefix_chromosomes(cis_frame['seqname'])
    cis_frame_r = pandas_to_dataframe(cis_frame)

    # Retrieve cis matrix from cimpl.
    cimpl_r = importr('cimpl')
    cis_matrix_r = cimpl_r.getCISMatrix(cimpl_obj, cis_frame_r)
    cis_matrix = dataframe_to_pandas(cis_matrix_r)

    # Extract scale information from cis matrix.
    scale_cols = [c for c in cis_matrix.columns if c.startswith('X')]
    cis_matrix_scales = cis_matrix[['id'] + scale_cols]

    # Melt matrix into long format.
    mapping = pd.melt(cis_matrix_scales, id_vars=['id'])
    mapping = mapping[['id', 'value']]
    mapping = mapping.rename(columns={'id': 'insertion_id', 'value': 'cis_id'})

    # Split cis_id column into individual entries (for entries
    # with multiple ids). Then drop any empty rows, as these
    # entries are empty cells in the matrix.
    mapping = _expand_column(mapping, col='cis_id', delimiter='|')
    mapping = mapping.ix[mapping['cis_id'] != '']

    return mapping


def merge_cis(cis_frame):
    """ Merge cis sites that are in fact the same, but appear multiple times
        with different peak_height locations.
    :param cis_frame:
    :return:
    """

    cols = ['chromosome', 'start', 'end', 'width', 'n_insertions', 'scale']
    return pd.DataFrame((grp.ix[grp['peak_height'].argmax()]
                         for _, grp in cis_frame.groupby(cols)))


def _expand_column(frame, col, delimiter):
    exp = pd.concat((_expand_row(row, col=col, delimiter=delimiter)
                     for _, row in frame.iterrows()), ignore_index=True)
    return exp[frame.columns]


def _expand_row(row, col, delimiter):
    row_dict = dict(row)

    if type(row[col]) == str:
        col_split = row[col].split(delimiter)
        row_dict[col] = col_split
    else:
        row_dict[col] = [row[col]]

    return pd.DataFrame(row_dict)
