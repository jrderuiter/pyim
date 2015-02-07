__author__ = 'Julian'


import readline  # Work-around for: ../lib/libreadline.so.6: undefined symbol: PC

import pandas as pd

from rpy2 import robjects
from rpy2.robjects.packages import importr

from .rpy import pandas_to_dataframe, dataframe_to_pandas

cimpl_r = importr('cimpl')


R_GENOMES = {
    'mm10': 'BSgenome.Mmusculus.UCSC.mm10'
}


def cimpl(insertions, scales, genome, system=None, specificity_pattern=None,
          n_iterations=100, chromosomes=None, verbose=False, threads=1):
    # Fill in chromosomes from data if not specified.
    if chromosomes is None:
        chromosomes = list(insertions['seqname'].unique())

    # Determine if system or specific pattern was specified.
    if specificity_pattern is not None:
        extra_args = {'specificity_pattern': specificity_pattern}
    elif system is not None:
        extra_args = {'system': system}
    else:
        raise ValueError('Either system or specificity pattern should be specified.')

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
    result = cimpl_r.doCimplAnalysis(
        _cimpl_frame(insertions), scales=scales, n_iterations=n_iterations, BSgenome=genome_obj,
        chromosomes=chromosomes, verbose=verbose, threads=threads, **extra_args)

    return result


def _cimpl_frame(insertions):
    # Extract and rename required columns.
    cimpl_frame = insertions.ix[:, ['name', 'seqname', 'location', 'sample']]
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


def cis(cimpl_obj, alpha=0.05, mul_test=True):
    cis_obj = cimpl_r.getCISs(cimpl_obj, alpha=alpha, mul_test=mul_test)

    # Convert cis to pandas and rename index.
    cis_frame = dataframe_to_pandas(cis_obj).reset_index()
    cis_frame.rename(columns={'index': 'id'}, inplace=True)

    # Convert columns to int types.
    for col in ['peak_location', 'start', 'end', 'width', 'n_insertions']:
        cis_frame[col] = cis_frame[col].astype(int)

    # Remove chr prefix from chromosomes.
    cis_frame['chromosome'] = cis_frame['chromosome'].str.replace('chr', '')

    return cis_frame


def cis_mapping(cimpl_obj, cis_frame):
    # Add cis_id as index to cis frame before passing to R,
    # ensures CIMPL uses cis id's instead of row indices.
    cis_r = cis_frame.copy()
    cis_r.index = cis_frame['id']
    cis_r['chromosome'] = _prefix_chromosomes(cis_r['chromosome'])
    cis_r = pandas_to_dataframe(cis_r)

    # Retrieve cis matrix from cimpl.
    cis_matrix_r = cimpl_r.getCISMatrix(cimpl_obj, cis_r)
    cis_matrix = dataframe_to_pandas(cis_matrix_r)

    # Extract scale information from cis matrix.
    scale_cols = [c for c in cis_matrix.columns if c.startswith('X')]
    cis_matrix_scales = cis_matrix[['id'] + scale_cols]

    # Melt matrix into long format.
    mapping = pd.melt(cis_matrix_scales, id_vars=['id'])[['id', 'value']]
    mapping = mapping.rename(columns={'id': 'insertion_id', 'value': 'cis_id'})

    # Split cis_id column into individual entries (for entries with multiple ids).
    # Then drop any empty rows, as these entries are empty cells in the matrix.
    mapping = _expand_column(mapping, col='cis_id', delim='|')
    mapping = mapping.ix[mapping['cis_id'] != '']

    return mapping


def _expand_column(frame, col, delim):
    exp = pd.concat((_expand_row(row, col=col, delim=delim)
                     for _, row in frame.iterrows()), ignore_index=True)
    return exp[frame.columns]


def _expand_row(row, col, delim):
    row_dict = dict(row)

    if type(row[col]) == str:
        col_split = row[col].split(delim)
        row_dict[col] = col_split
    else:
        row_dict[col] = [row[col]]

    return pd.DataFrame(row_dict)
