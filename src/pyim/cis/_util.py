import toolz

import numpy as np
import pandas as pd


def annotate_cis_strand(cis, insertions, min_homogeneity):
     # Determine strand of cis sites.
     func = toolz.curry(_cis_strand, min_homogeneity=min_homogeneity)
     cis_strand = insertions.groupby('cis_id').apply(func)

     # Merge with cis annotation
     cis = pd.merge(cis, cis_strand.reset_index(), on='cis_id')

     return cis


def _cis_strand(insertions, min_homogeneity):
     strand_mean = insertions.strand.mean()
     strand = int(np.sign(strand_mean))

     if strand != 0:
         homogeneity = (insertions.strand == strand).sum() / len(insertions)
     else:
         homogeneity = 0.5

     if homogeneity < min_homogeneity:
         strand = 0

     return pd.Series(dict(strand=strand,
                           strand_mean=strand_mean,
                           strand_homogeneity=homogeneity))