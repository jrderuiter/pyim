import pandas as pd
import numpy as np

from pyim.model import InsertionSet


def filter_matches(matches,
                   closest=False,
                   blacklist=None,
                   blacklist_col='gene_id'):
    """Filter matches for closest gene and/or blacklisted genes."""

    if closest:
        grouped = matches.groupby(['id'])
        mask = grouped['gene_distance'].transform(lambda x: x == x.min())

        matches = matches.loc[mask]

    if blacklist is not None:
        matches = matches.loc[~matches[blacklist_col].isin(set(blacklist))]

    return matches


def annotate_matches(matches, insertions, genes):
    """Annotate matches with gene distance/orientation."""

    # TODO: Handle different chromosome case in distance?
    #   (Even if this not expected in the input.)

    # Add gene/insertion info to matches.
    gene_col_mapping = {
        'gene_id': 'gene_id',
        'gene_name': 'gene_name',
        'start': 'gene_start',
        'end': 'gene_end',
        'strand': 'gene_strand'
    }

    gene_info = (genes[list(gene_col_mapping.keys())]
                 .rename(columns=gene_col_mapping))  # yapf: disable

    gene_info['gene_strand'] = gene_info['gene_strand'].map({'+': 1, '-': -1})

    ins_col_mapping = {
        'id': 'id',
        'position': 'ins_position',
        'strand': 'ins_strand'
    }

    ins_info = (insertions.values[list(ins_col_mapping.keys())]
                .rename(columns=ins_col_mapping))  # yapf: disable

    merged = pd.merge(matches, gene_info, on='gene_id', how='left')
    merged = pd.merge(merged, ins_info.drop_duplicates(), on='id', how='left')

    # Calculate distances.
    start_dist = (merged['gene_start'] - merged['ins_position']).abs()
    end_dist = (merged['gene_end'] - merged['ins_position']).abs()

    merged['gene_distance'] = np.minimum(start_dist, end_dist)

    within_mask = ((merged['gene_start'] < merged['ins_position']) &
                   (merged['ins_position'] < merged['gene_end']))
    merged.loc[within_mask, 'gene_distance'] = 0

    # Determine (relative) orientation.
    same_orientation = merged['gene_strand'] == merged['ins_strand']
    merged['gene_orientation'] = same_orientation.map({True: 1, False: -1})

    return merged.reindex(columns=list(matches.columns) +
                          ['gene_name', 'gene_distance', 'gene_orientation'])


# def annotate_matched_insertions(insertions, genes):
#     """Annotate matched insertions with gene distance/orientation."""

#     annotated = annotate_matches(
#         insertions.values[['id', 'gene_id']].drop_duplicates(),
#         insertions=insertions,
#         genes=genes)

#     # Drop existing annotation columns.
#     annot_cols = ['gene_name', 'gene_distance', 'gene_orientation']
#     ins_values = insertions.values.drop(annot_cols, errors='ignore', axis=1)

#     merged = pd.merge(ins_values, annotated, on=['id', 'gene_id'], how='left')

#     return InsertionSet(merged)
