import toolz
import numpy as np
import pandas as pd
import scipy.cluster.hierarchy as sch
import scipy.spatial.distance as ssd


def merge_within_distance(insertions, max_dist=2000, agg_funcs=None):
    clustered = cluster_insertions(insertions, max_dist=max_dist)
    return merge_insertions(clustered, by='cluster', agg_funcs=agg_funcs)


def cluster_insertions(insertions, max_dist=2000, method='single'):
    prev_n_clusters = 0

    clustered_grps = []
    for _, group in insertions.groupby(['chrom', 'barcode', 'strand']):
        clusters = _cluster_group(group, max_dist, method)
        clusters += prev_n_clusters

        group = group.copy()
        group['cluster'] = clusters
        clustered_grps.append(group)

        prev_n_clusters = np.max(clusters)

    return pd.concat(clustered_grps, ignore_index=True)


def _cluster_group(insertions, max_dist, method):
    if len(insertions) == 1:
        clusters = np.array([1], dtype=np.int32)
    else:
        dists = genomic_distance(insertions)
        z = sch.linkage(dists, method=method)
        clusters = sch.fcluster(z, criterion='distance', t=max_dist)

    return clusters


def genomic_distance(insertions):
    # Sanity check insertions (for debugging).
    assert(insertions['chrom'].nunique() == 1)
    assert(insertions['barcode'].nunique() == 1)
    assert(insertions['strand'].nunique() == 1)

    # Calculate 1d distances.
    loc = insertions['position']
    loc_2d = np.vstack([loc, np.zeros_like(loc)]).T
    dist = ssd.pdist(loc_2d, lambda u, v: np.abs(u-v).sum())

    return dist


def merge_insertions(insertions, by='cluster', agg_funcs=None):
    # TODO: use weighted median.
    default_agg = {'id': 'first', 'chrom': 'first', 'position': 'median',
                   'strand': 'first', 'barcode': 'first'}
    agg_funcs = toolz.merge(default_agg, agg_funcs or {})

    col_order = [c for c in insertions.columns if c in agg_funcs]
    merged = insertions.groupby(by).agg(agg_funcs)[col_order]

    return merged
