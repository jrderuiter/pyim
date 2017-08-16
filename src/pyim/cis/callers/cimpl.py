import pandas as pd

import readline
from rpy2 import robjects
from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import DataFrame as RDataFrame, StrVector, IntVector
from pyim.util.rpy2 import pandas_to_dataframe, dataframe_to_pandas

from pyim.model import Insertion, CisSite
from pyim.util import add_prefix, remove_prefix

from .base import CisCaller, CisCallerCommand
from ..util import assign_strand, invert_otm_mapping

R_GENOMES = {'mm10': 'BSgenome.Mmusculus.UCSC.mm10'}


class CimplCisCaller(CisCaller):
    def __init__(self,
                 genome='mm10',
                 scales=(10000, 30000),
                 chromosomes=None,
                 alpha=0.05,
                 pattern=None,
                 lhc_method='none',
                 iterations=1000,
                 threads=1,
                 min_strand_homogeneity=0.75):
        super().__init__()

        # Default to numbered mouse chromosomes + X.
        if chromosomes is None:
            chromosomes = [str(i) for i in range(1, 20)] + ['X']

        # Add 'chr' prefix to chromosomes if missing.
        chromosomes = add_prefix(chromosomes, prefix='chr')

        self._genome = genome
        self._scales = scales
        self._chromosomes = chromosomes
        self._alpha = alpha
        self._pattern = pattern
        self._lhc_method = lhc_method
        self._iterations = iterations
        self._threads = threads
        self._min_strand_homogeneity = min_strand_homogeneity

        self.__cimpl = None

    @property
    def _cimpl(self):
        if self.__cimpl is None:
            self.__cimpl = importr('cimpl')
        return self.__cimpl

    def call(self, insertions):
        """Runs CIMPL on insertions."""

        # Convert insertions to cimpl frame.
        ins_frame = self._insertions_to_cimpl(insertions)

        # Load genome object from R.
        genome_obj = self._load_genome(self._genome)

        # Check if contig_depth is present (if doing hop exclusion).
        if self._lhc_method == 'exclude' and 'contig_depth' not in ins_frame:
            raise ValueError('Insertion depth is needed for lhc exclusion')

        # Run CIMPL!
        cimpl_result = self._cimpl.doCimplAnalysis(
            pandas_to_dataframe(ins_frame),
            scales=robjects.vectors.IntVector(self._scales),
            n_iterations=self._iterations,
            lhc_method=self._lhc_method,
            cores=self._threads,
            BSgenome=genome_obj,
            chromosomes=robjects.vectors.StrVector(self._chromosomes),
            verbose=1)

        # Extract cis sites and mapping.
        cis_sites = self._extract_cis(cimpl_result, alpha=self._alpha)
        cis_mapping = self._extract_mapping(cimpl_result, cis_sites)

        # Determine strandedness of cis_sites using inseriions.
        ins_mapping = invert_otm_mapping(cis_mapping)
        cis_sites = list(
            assign_strand(
                cis_sites,
                insertions,
                ins_mapping,
                min_homogeneity=self._min_strand_homogeneity))

        return cis_sites, cis_mapping

    def _insertions_to_cimpl(self, insertions):
        # Convert insertions to frame representation.
        ins_frame = Insertion.to_frame(insertions)

        # Extract and rename required columns.
        column_map = {
            'id': 'id',
            'chromosome': 'chr',
            'position': 'location',
            'sample': 'sampleID'
        }

        cimpl_ins = (ins_frame[list(column_map.keys())]
                     .rename(columns=column_map))

        # Add chr prefix.
        cimpl_ins['chr'] = add_prefix(cimpl_ins['chr'], prefix='chr')

        # Add depth if present.
        if 'support' in ins_frame:
            cimpl_ins['contig_depth'] = ins_frame['support']

        elif 'depth_unique' in ins_frame:
            cimpl_ins['contig_depth'] = ins_frame['depth_unique']

        return cimpl_ins

    def _load_genome(self, genome):
        # Lookup R package for genome.
        try:
            genome_pkg = R_GENOMES[genome]
        except KeyError:
            raise ValueError('Unsupported genome {}'.format(genome))

        # Import package and extract genome object.
        bs_genome = importr(genome_pkg)
        genome_obj = bs_genome.Mmusculus

        return genome_obj

    def _extract_cis(self, cimpl_obj, alpha=0.05):
        # Extract CIS frame from R into a pandas dataframe.
        cis_obj = self._cimpl.getCISs(cimpl_obj, alpha=alpha, mul_test=True)

        cis_frame = dataframe_to_pandas(cis_obj).reset_index()
        cis_frame.rename(
            columns={'index': 'cis_id',
                     'chromosome': 'seqname'}, inplace=True)

        # Clean-up converted dataframe (datatypes, prefixes + column names).
        for col in ['peak_location', 'start', 'end', 'width', 'n_insertions']:
            cis_frame[col] = cis_frame[col].astype(int)

        cis_frame['seqname'] = remove_prefix(
            cis_frame['seqname'], prefix='chr')

        cis_frame = cis_frame.rename(columns={
            'cis_id': 'id',
            'seqname': 'chromosome',
            'peak_location': 'position',
            'peak_height': 'height',
            'p_value': 'pvalue'
        })

        # Merge cis sites that are in fact the same, but appear multiple
        # times with different height locations.
        cols = ['chromosome', 'start', 'end', 'width', 'n_insertions', 'scale']
        cis_frame = pd.DataFrame((grp.ix[grp['height'].argmax()]
                                  for _, grp in cis_frame.groupby(cols)))

        # For now, set strand to None.
        cis_frame['strand'] = None

        # Convert to CisSite objects using a subset of the columns.
        cis_frame_subset = cis_frame[[
            'id', 'chromosome', 'position', 'start', 'end', 'scale', 'pvalue',
            'n_insertions', 'height', 'width', 'strand'
        ]]
        cis_sites = list(CisSite.from_frame(cis_frame_subset))

        return cis_sites

    def _extract_mapping(self, cimpl_obj, cis_sites):
        # Convert CIS sites to frame format.
        cis_frame = CisSite.to_frame(cis_sites)

        # Convert to R representation for cimpl.
        chr_with_prefix = add_prefix(cis_frame['chromosome'], prefix='chr')

        r_base = importr('base')
        cis_frame_r = RDataFrame({
            'id': r_base.I(StrVector(cis_frame['id'])),
            'chromosome': r_base.I(StrVector(chr_with_prefix)),
            'scale': StrVector(cis_frame['scale']),
            'start': IntVector(cis_frame['start']),
            'end': IntVector(cis_frame['end'])
        })
        cis_frame_r.rownames = StrVector(cis_frame['id'])

        # Retrieve cis matrix from cimpl.
        cis_matrix_r = self._cimpl.getCISMatrix(cimpl_obj, cis_frame_r)
        cis_matrix = dataframe_to_pandas(cis_matrix_r)

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

        mapping_dict = {
            ins_id: set(grp['cis_id'])
            for ins_id, grp in mapping.groupby('insertion_id')
        }

        return mapping_dict


class CimplCisCallerCommand(CisCallerCommand):

    name = 'cimpl'

    def configure(self, parser):
        super().configure(parser)

        parser.add_argument('--pattern', required=True)

        parser.add_argument('--genome', default='mm10')
        parser.add_argument('--scales', default=(10000, 30000),
                            nargs='+', type=int)  # yapf: disable
        parser.add_argument('--chromosomes', default=None, nargs='+')
        parser.add_argument('--alpha', default=0.05, type=float)
        parser.add_argument('--lhc_method', default='exclude')
        parser.add_argument('--iterations', default=1000, type=int)
        parser.add_argument('--threads', default=1, type=int)
        parser.add_argument(
            '--min_strand_homogeneity', default=0.75, type=float)

    def run(self, args):
        # Read insertions.
        insertions = list(Insertion.from_csv(args.insertions, sep='\t'))

        # Call CIS sites.
        caller = CimplCisCaller(
            pattern=args.pattern,
            genome=args.genome,
            chromosomes=args.chromosomes,
            alpha=args.alpha,
            lhc_method=args.lhc_method,
            iterations=args.iterations,
            threads=args.threads,
            min_strand_homogeneity=args.min_strand_homogeneity)

        cis_sites, cis_mapping = caller.call(insertions=insertions)

        # Annotate insertions and write outputs.
        annotated_ins = self._annotate_insertions(insertions, cis_mapping)
        self._write_outputs(annotated_ins, cis_sites, args)


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
