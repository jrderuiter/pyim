import pandas as pd

import readline
from rpy2 import robjects
from rpy2.robjects.packages import importr

from pyim.model import InsertionSet, CisSet
from pyim.util.rpy2 import pandas_to_dataframe, dataframe_to_pandas

from .base import CisCaller

R_GENOMES = {'mm10': 'BSgenome.Mmusculus.UCSC.mm10'}

_cimpl = None


def _get_cimpl():
    global _cimpl

    if _cimpl is None:
        _cimpl = importr('cimpl')

    return _cimpl


class CimplCisCaller(CisCaller):
    def __init__(self,
                 genome='mm10',
                 scales=(10000, 30000),
                 chromosomes=None,
                 alpha=0.05,
                 pattern=None,
                 lhc_method='none',
                 iterations=1000,
                 threads=1):
        super().__init__()

        # Default to numbered mouse chromosomes + X.
        if chromosomes is None:
            chromosomes = [str(i) for i in range(1, 20)] + ['X']

        self._genome = genome
        self._scales = scales
        self._chromosomes = chromosomes
        self._alpha = alpha
        self._pattern = pattern
        self._lhc_method = lhc_method
        self._iterations = iterations
        self._threads = threads

    def call(self, insertions):
        """Runs CIMPL on insertions."""

        # Convert insertions to cimpl format.
        cimpl_ins = self._insertions_to_cimpl(insertions)
        chromosomes = _add_prefix(self._chromosomes, prefix='chr')

        # Check if contig_depth is present (if doing hop exclusion).
        if self._lhc_method == 'exclude' and 'contig_depth' not in cimpl_ins:
            raise ValueError('Insertion depth is needed for lhc exclusion')

        # Run CIMPL!
        cimpl = _get_cimpl()

        result = CimplResult(
            cimpl.doCimplAnalysis(
                pandas_to_dataframe(cimpl_ins),
                BSgenome=self._load_genome(self._genome),
                chromosomes=robjects.vectors.StrVector(chromosomes),
                scales=robjects.vectors.IntVector(self._scales),
                n_iterations=self._iterations,
                lhc_method=self._lhc_method,
                cores=self._threads,
                verbose=1))

        # Extract cis sites and mapping.
        cis_sites = result.get_cis_sites()
        cis_mapping = result.get_cis_mapping()

        return cis_sites, cis_mapping

    def _insertions_to_cimpl(self, insertions):
        """Converts InsertionSet to CIMPL frame representation."""

        # Extract and rename required columns.
        column_map = {
            'id': 'id',
            'chromosome': 'chr',
            'position': 'location',
            'sample': 'sampleID',
            'support': 'contig_depth'
        }

        cimpl_ins = (insertions.values.rename(columns=column_map)[[
            'id', 'chr', 'location', 'sampleID', 'contig_depth'
        ]])

        # Add chr prefix.
        cimpl_ins['chr'] = _add_prefix(cimpl_ins['chr'], prefix='chr')

        return cimpl_ins

    def _load_genome(self, genome):
        """Loads given reference genome."""

        # Lookup R package for genome.
        try:
            genome_pkg = R_GENOMES[genome]
        except KeyError:
            raise ValueError('Unsupported genome {}'.format(genome))

        # Import package and extract genome object.
        bs_genome = importr(genome_pkg)
        genome_obj = bs_genome.Mmusculus

        return genome_obj


class CimplResult(object):
    """Wrapper class for CIMPL R results."""

    def __init__(self, result):
        self._result = result

    def get_cis_sites(self, alpha=0.05):
        """Extracts CIS sites from the CIMPL result."""

        # Get CIS object from CIMPL.
        cimpl = _get_cimpl()
        cis_r = cimpl.getCISs(self._result, alpha=alpha, mul_test=True)

        # Extract CIS data in DF format.
        cis_data = (dataframe_to_pandas(cis_r)
                    .reset_index()
                    .rename(columns={'index': 'id'}))  # yapf: disable

        # Map and extract required columns.
        col_mapping = {
            'cis_id': 'id',
            'seqname': 'chromosome',
            'peak_location': 'position',
            'peak_height': 'height',
        }

        cols = [
            'id', 'chromosome', 'position', 'start', 'end', 'scale', 'p_value',
            'n_insertions', 'height', 'width', 'strand'
        ]
        cis_data = cis_data.rename(columns=col_mapping).reindex(columns=cols)

        # Fix column datatypes.
        for col in ['position', 'start', 'end', 'width', 'n_insertions']:
            cis_data[col] = cis_data[col].astype(int)

        # For now, set strand to None.
        cis_data['strand'] = None

        return cis_data

    def get_cis_mapping(self, alpha=0.05):
        """Extracts insertion-to-cis mapping."""

        # Get CIS matrix object from CIMPL.
        cimpl = _get_cimpl()

        cis_r = cimpl.getCISs(self._result, alpha=alpha, mul_test=True)
        cis_matrix = dataframe_to_pandas(
            cimpl.getCISMatrix(self._result, cis_r))

        # Extract scale information from cis matrix.
        scale_cols = [c for c in cis_matrix.columns if c.startswith('X')]
        cis_matrix_scales = cis_matrix[['id'] + scale_cols]

        # Melt matrix into long format.
        mapping = pd.melt(cis_matrix_scales, id_vars=['id'])
        mapping = mapping[['id', 'value']]
        mapping = mapping.rename(
            columns={'id': 'insertion_id',
                     'value': 'cis_id'})

        # Split cis_id column into individual entries (for entries
        # with multiple ids). Then drop any empty rows, as these
        # entries are empty cells in the matrix.
        mapping = mapping.ix[mapping['cis_id'] != '']
        mapping = expand_column(mapping, col='cis_id', delimiter='|')

        return mapping


def expand_column(frame, col, delimiter):
    exp = pd.concat(
        (_expand_row(row, col=col, delimiter=delimiter)
         for _, row in frame.iterrows()),
        ignore_index=True)  # yapf: disable
    return exp[frame.columns]


def _expand_row(row, col, delimiter):
    row_dict = dict(row)

    if type(row[col]) == str:
        col_split = row[col].split(delimiter)
        row_dict[col] = col_split
    else:
        row_dict[col] = [row[col]]

    return pd.DataFrame(row_dict)


def _add_prefix(values, prefix='chr'):
    """Adds prefix to (str) values."""

    def _prefix(value):
        if not value.startswith(prefix):
            return prefix + value
        return value

    return [_prefix(value) for value in values]


def _remove_prefix(values, prefix='chr'):
    """Removes prefix from (str) values."""

    def _remove(value):
        if value.startswith(prefix):
            return value[len(prefix):]
        return value

    return [_prefix(value) for value in values]
