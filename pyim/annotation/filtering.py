
def select_closest(insertions, id_col='id', dist_col='distance'):
    """Selects genes that are closest to the annotated insertions.

    Args:
        insertions (pandas.DataFrame): Annotated insertions that are to
            be filtered. The frame is expected to contain at least the
            following columns: id, position, strand, *dist_col*.
        id_col (str): Name of the column containing the id of the insertion.
        dist_col (str): Name of the column containing the distance to
            the gene or feature. Can be added using the add_metadata function.

    Returns:
        pandas.DataFrame: Filtered annotated insertions, which have been
            reduced to only include the genes closest to the insertions.

    """

    def _is_closest(x):
        abs_dist = x[dist_col].abs()
        return x.ix[abs_dist == abs_dist.min()]

    return (insertions.groupby(id_col)
            .apply(_is_closest)
            .reset_index(drop=True))


def filter_blacklist(insertions, blacklist, gene_col='gene_name'):
    """Filters annotations that assign insertions to blacklisted genes.

    Args:
        insertions (pandas.DataFrame): Annotated insertions that are to
            be filtered. The frame is expected to contain at least the
            following columns: id, position, strand, *gene_id_col*.
        blacklist (list[str]): List of blacklisted gene ids to filter.
        gene_col (str): Name of the column containing the id of the genes.

    Returns:
        pandas.DataFrame: Filtered annotated insertions, which have been
            reduced remove blacklisted genes.

    """

    mask = insertions[gene_col].isin(set(blacklist))
    return insertions.ix[~mask]
