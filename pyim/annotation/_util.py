
def select_closest(frame, id_col='id', col='distance'):
    def _is_closest(x):
        abs_dist = x[col].abs()
        return x.ix[abs_dist == abs_dist.min()]

    return (frame.groupby(id_col)
            .apply(_is_closest)
            .reset_index(drop=True))


def feature_distance(feature, location, stranded=True):
    start, end = feature['start'], feature['end']

    if start <= location <= end:
        dist = 0
    elif location > end:
        dist = location - end
    else:
        dist = location - start

    if stranded:
        dist *= numeric_strand(feature['strand'])

    return dist


def numeric_strand(strand):
    """Convert strand to numeric representation."""

    return 1 if strand == '+' else -1


