
import pandas as pd, numpy as np

columnNames = ['match', 'mismatch', 'repmatch', 'Ns',
               'qGapCount', 'qGapBases',
               'tGapCount', 'tGapBases',
               'strand',
               'qName', 'qSize', 'qStart', 'qEnd',
               'tName', 'tSize', 'tStart', 'tEnd',
               'blockCount', 'blockSizes',
               'qStarts', 'tStarts']


def read_psl(filePath, onlyPutative=False):
    pslReads = pd.read_table(filePath, header=None, skiprows=5, names=columnNames, index_col='qName')
    if onlyPutative:
        pslReads = select_putative_reads(pslReads)
    return pslReads


def select_putative_reads(pslReads):
    def selmax(group):
        maxInd = np.argmax(group['match'])
        return group.ix[maxInd]

    grp = pslReads.groupby(level=0)
    return grp.apply(selmax)