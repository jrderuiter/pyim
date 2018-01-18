import gzip
import pickle

import pytest

from pyim.cis.callers import cimpl


class CimplCisCaller(object):
    pass


class TestCimplResult(object):
    @pytest.fixture(scope='class')
    def result(self):
        result_path = pytest.helpers.data_path('sb/cimpl.pkl.gz')

        with gzip.open(str(result_path), 'rb') as file_:
            result_obj = pickle.load(file_)

        return cimpl.CimplResult(result_obj)

    def test_get_cis_sites(self, result):
        print(result.get_cis_sites())

    def test_get_cis_mapping(self, result):
        print(result.get_cis_mapping())

    def run(self, args):
        # Read insertions.
        insertions = InsertionSet.from_csv(args.insertions, sep='\t')

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
