
import numpy
import pandas
from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import single, fcluster


def cluster_insertions(insertions, max_dist=10):
    clustered_insertions = []

    for _, ins_grp in insertions.groupby(['barcode', 'chromosome', 'strand']):
        if len(ins_grp) == 1:
            clustered_grp = ins_grp
        else:
            dist = _insertion_distances(ins_grp)
            clust = fcluster(single(dist), t=max_dist, criterion='distance')
            clustered_grp = _apply_clustering(ins_grp, clust)
        clustered_insertions.append(clustered_grp)

    return pandas.concat(clustered_insertions, ignore_index=True)


def _insertion_distances(insertions):
    loc = insertions['location']

    loc_2d = numpy.vstack([loc, numpy.zeros_like(loc)]).T
    dist = pdist(loc_2d, lambda u,v: numpy.abs(u-v).sum())

    return dist


def _apply_clustering(insertions, clusters):
    clustered_insertions = []
    for i in range(1, max(clusters) + 1):
        cluster_ins = insertions.ix[clusters == i]

        # Create new putative insertion from cluster by selecting
        # the insertion with the highest ULP as the putative insertion.
        max_ulp_ind = cluster_ins['unique_lp'].argmax()

        putative_ins = cluster_ins.ix[max_ulp_ind].copy()
        putative_ins['lp'] = cluster_ins['lp'].sum()
        putative_ins['unique_lp'] = cluster_ins['unique_lp'].sum()

        clustered_insertions.append(putative_ins)

    return pandas.DataFrame(clustered_insertions)     # Return insertions as frame
