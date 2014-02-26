
import numpy, pandas
from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import single, fcluster


def cluster_insertions(insertions, dist_t=10):
    ins_clust = insertion_clustering(insertions)
    return putative_insertions(ins_clust)


def insertion_clustering(insertions, dist_t=10):
    clustered, clust_offset = [], 0

    for _, group in insertions.groupby(['chromosome', 'strand', 'barcode']):
        if len(group) > 1:
            locations_2d = numpy.vstack([group['location'], numpy.zeros_like(group['location'])]).T
            dist = pdist(locations_2d, lambda u,v: numpy.abs(u-v).sum())
            z = single(dist)
            clust = fcluster(z, t=dist_t, criterion='distance')

            clust_labels = clust + clust_offset
            clust_offset += clust.max()
        else:
            clust_labels = 1 + clust_offset
            clust_offset += 1

        group['clust'] = clust_labels
        clustered.append(group)

    return pandas.concat(clustered)


def putative_insertions(ins_clust):
    insertions = []

    for _, cluster in ins_clust.groupby('clust'):
        #if np.sum(cluster['unique_lp'] > 3) > 1:
        #    print cluster
        putative_ins = cluster.iloc[cluster['unique_lp'].argmax()]
        putative_ins['lp'] = cluster['lp'].sum()
        putative_ins['unique_lp'] = cluster['unique_lp'].sum()
        insertions.append(putative_ins)

    return pandas.DataFrame(insertions)