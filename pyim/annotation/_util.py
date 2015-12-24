
def get_closest(frame, id_col='insertion_id', distance_col='distance'):
    def _is_closest(x):
        abs_dist = x[distance_col].abs()
        return x.ix[abs_dist == abs_dist.min()]

    return (frame.groupby(id_col)
            .apply(_is_closest)
            .reset_index(drop=True))


def feature_distance(start, end, location):
    if start <= location <= end:
        return 0
    elif location > end:
        return location - end
    else:
        return location - start
