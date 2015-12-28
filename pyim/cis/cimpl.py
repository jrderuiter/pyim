import pandas as pd

import readline
from rpy2 import robjects
from rpy2.robjects.packages import importr

from pyim.util.rpy2 import pandas_to_dataframe, dataframe_to_pandas


R_GENOMES = {
    'mm10': 'BSgenome.Mmusculus.UCSC.mm10'
}


def map_insertions(insertions, scales, genome, alpha=0.05, **kwargs):
    """Maps given insertions to CISs using CIMPL."""

    # Convert insertion to cimpl format.
    cimpl_ins = convert_to_cimpl(insertions)

    # Run cimpl.
    cimpl_result = cimpl(cimpl_ins, scales, genome, **kwargs)

    # Extract cis sites and mapping.
    cis = extract_cis(cimpl_result, alpha=alpha)
    mapping = extract_mapping(cimpl_result, cis)

    return cis, mapping


def cimpl(insertions, scales, genome, system=None, pattern=None,
          lhc_method='none', iterations=1000, chromosomes=None,
          verbose=False, threads=1):
    """Runs CIMPL on insertions (in CIMPL format)."""

    # Fill in chromosomes from data if not specified.
    if chromosomes is None:
        chromosomes = list(insertions['chr'].unique())

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

    # Convert scales to IntVector if supplied as list.
    if type(scales) == list:
        scales = robjects.vectors.IntVector(scales)

    # Load genome object from R.
    genome_obj = _load_genome(genome)

    # Check if contig_depth is present (if doing hop exclusion).
    if lhc_method == 'exclude' and 'contig_depth' not in insertions:
        raise ValueError('Insertion depth is needed for lhc exclusion')

    # Run CIMPL!
    cimpl_r = importr('cimpl')
    cimpl_obj = cimpl_r.doCimplAnalysis(
        pandas_to_dataframe(insertions),
        scales=scales, n_iterations=iterations,
        lhc_method=lhc_method, threads=threads, BSgenome=genome_obj,
        chromosomes=robjects.vectors.StrVector(chromosomes),
        verbose=verbose, **extra_args)

    return cimpl_obj


def convert_to_cimpl(insertions):
    # Extract and rename required columns.
    cimpl_ins = insertions.ix[:, ['id', 'chrom', 'position', 'sample']]
    cimpl_ins.columns = ['id', 'chr', 'location', 'sampleID']

    if 'depth_unique' in insertions:
        cimpl_ins['contig_depth'] = insertions['depth_unique']

    # Add 'chr' prefix to the chromosome names if needed.
    cimpl_ins['chr'] = _prefix_chromosomes(cimpl_ins['chr'])

    return cimpl_ins


def _prefix_chromosomes(series, prefix='chr'):
    # Add 'chr' prefix to the chromosome names if needed.
    if len(series) > 0 and not series.iloc[0].startswith('chr'):
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


def extract_cis(cimpl_obj, alpha=0.05, mul_test=True):
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

    # Rename and reshuffle cis columns.
    cis_frame = cis_frame.rename(
        columns={'seqname': 'chrom',
                 'peak_location': 'position',
                 'peak_height': 'height'})

    cis_frame = cis_frame[['cis_id', 'chrom', 'position', 'scale',
                           'n_insertions', 'p_value', 'start', 'end',
                           'height', 'width']]

    return cis_frame


def extract_mapping(cimpl_obj, cis_frame):
    # Add cis_id as index to cis frame before passing to R,
    # ensures CIMPL uses cis id's instead of row indices.
    cis_frame = cis_frame.copy()
    cis_frame.set_index('cis_id', drop=False, inplace=True)
    cis_frame['chromosomes'] = _prefix_chromosomes(cis_frame['chrom'])
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


# def cis_strandedness(insertions, min_homogeneity):
#     strand_mean = insertions.strand.mean()
#     strand = int(np.sign(strand_mean))
#
#     if strand != 0:
#         homogeneity = (insertions.strand == strand).sum() / len(insertions)
#     else:
#         homogeneity = 0.5
#
#     if homogeneity < min_homogeneity:
#         strand = 0
#
#     return pd.Series(dict(strand=strand,
#                           strand_mean=strand_mean,
#                           strand_homogeneity=homogeneity))
#
#     # Determine strand of cis sites.
#     strand_func = curry(_strandedness, min_homogeneity=args.strand_homogeneity)
#     cis_strand = insertions.groupby('cis_id').apply(strand_func)
#     cis = pd.merge(cis, cis_strand.reset_index(), on='cis_id')
