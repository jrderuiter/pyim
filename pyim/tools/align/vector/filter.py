__author__ = 'Julian'


def identity_filter(aln, min_identity):
    return aln.identity >= min_identity


def coverage_filter(aln, min_coverage):
    return aln.coverage >= min_coverage
