from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
from builtins import (bytes, dict, int, list, object, range, str,
                      ascii, chr, hex, input, next, oct, open,
                      pow, round, super, filter, map, zip)

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster, complete, average

LINKAGE_MAP = {
    'complete': complete,
    'average': average
}


def cluster_frame(frame, dist_func, groupby=None,
                  linkage='complete', criterion='distance', t=0):
    # Lookup linkage type.
    try:
        linkage = LINKAGE_MAP[linkage]
    except KeyError:
        raise ValueError('Unknown linkage type {}'.format(linkage))

    # Group by columns if any are given.
    if groupby is None:
        groups = [(None, frame)]
    else:
        # Replace NaNs to avoid dropping entries.
        frame = frame.fillna('NaN')
        groups = frame.groupby(groupby)

    # Determine clusters and use to sub-group frame.
    for _, grp in groups:
        if len(grp) == 1:
            yield grp.replace('NaN', np.nan)
        else:
            dists = dist_func(grp)
            clusters = fcluster(complete(dists), t=t, criterion=criterion)
            for _, cluster_grp in grp.groupby(clusters, group_keys=False):
                if groupby is not None:
                    cluster_grp = cluster_grp.replace('NaN', np.nan)
                yield cluster_grp


def cluster_frame_merged(frame, dist_func, merge_func, **kwargs):
    # Otherwise we merge each group into a single row (series)
    # and return the summarized dataframe.
    groups = list(cluster_frame(frame, dist_func, **kwargs))

    frame = pd.DataFrame.from_records(
        (merge_func(grp) for grp in groups),
        columns=frame.columns)

    return frame


