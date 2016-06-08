import itertools
import re

import toolz
import pandas as pd
from intervaltree import IntervalTree
from scipy.stats import poisson
from statsmodels.stats.multitest import multipletests


def build_trees(insertions):
    trees = {}

    for chrom, grp in insertions.groupby('chrom'):
        intervals = zip(grp['position'], grp['position'] + 1, grp['id'])
        trees[chrom] = IntervalTree.from_tuples(intervals)

    return trees


def count_pattern(record, pattern=None):
    regex = re.compile(pattern)
    return sum((1 for match in regex.finditer(record.seq)))


def count_matches(seq, regex):
    return sum((1 for match in regex.finditer(seq)))


def generate_windows(insertions, window_size):
    half_size = window_size // 2

    # Generate list of windows for all insertions.
    windows = (zip((chrom for _ in range(len(grp))),
                   grp['position'] - half_size,
                   grp['position'] + half_size)
               for chrom, grp in insertions.groupby('chrom'))
    windows = itertools.chain.from_iterable(windows)

    # Yield from windows.
    for window in windows:
        yield window


def calc_significance(insertions, reference, window_size,
                      pattern=None, chromosomes=None, total=None):
    if chromosomes is None:
        chromosomes = reference.keys()

    if pattern is not None:
        regex = re.compile(pattern)
        func = toolz.curry(count_matches, regex=regex)
    else:
        func = len

    if total is None:
        total = sum((func(reference[c][0:len(reference[c])].seq)
                     for c in chromosomes))

    # Subset insertions to chromosomes:
    insertions = insertions.ix[
        insertions['chrom'].isin(chromosomes)]

    # Build lookup trees for insertions.
    trees = build_trees(insertions)

    def _calc_for_window(window):
        chrom, start, end = window

        # Calculate occurrence for region.
        n_region = func(reference[chrom][int(start):int(end)].seq)

        # Calculate p-value.
        x = len(trees[chrom][start:end])
        mu = len(insertions) * (n_region / total)

        p_val = poisson.sf(x, mu=mu, loc=1)

        return chrom, start, end, p_val

    # Generate windows.
    windows = generate_windows(insertions, window_size=window_size)

    # Generate result.
    res = pd.DataFrame.from_records(
        (_calc_for_window(w) for w in windows),
        columns=['chrom', 'start', 'end', 'p_val'])
    res['p_val_corr'] = multipletests(res['p_val'], method='bonferroni')[1]

    return res


# result = calc_significance(insertions, ref, window_size=10000,
#                            pattern='(AT|TA)', chromosomes=chroms,
#                            total=genome_ta)
# result.query('p_val_corr < 0.05')
