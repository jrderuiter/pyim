
import pandas as pd, numpy as np

columnNames = ['match', 'mismatch', 'repmatch', 'num_n',
               'q_gap_count', 'q_gap_bases',
               't_gap_count', 't_gap_bases',
               'strand',
               'q_name', 'q_size', 'q_start', 'q_end',
               't_name', 't_size', 't_start', 't_end',
               'block_count', 'block_sizes',
               'q_starts', 't_starts']


def read_psl(filePath, onlyPutative=False):
    hits = pd.read_table(filePath, header=None, skiprows=5, names=columnNames, index_col='q_name')
    if onlyPutative:
        hits = select_putative_reads(hits)
    return hits


def select_putative_reads(hits):
    def selmax(group):
        maxInd = np.argmax(group['match'])
        return group.ix[maxInd]

    grp = hits.groupby(level=0)
    return grp.apply(selmax)