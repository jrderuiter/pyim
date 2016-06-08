import pandas as pd
from pyim.util.tabix import GtfFile, GtfFrame

from .util import numeric_strand


def add_metadata(insertions, gtf):
    """Adds metadata to annotated insertions.

    Adds extra metadata to already annotated insertions. This metadata
    currently includes the following information: distance to the gene
    ('distance' column) and relative orientation ('orientation' column).

    Args:
        insertions (pandas.DataFrame): Annotated insertions for which metadata
            should be added. The frame is expected to contain at least the
            following columns: id, position, strand, gene_id.
        gtf (str or GtfFile):  Path to gtf file containing gene features.
            Alternatively, a GtfFile object may also be given instead of a path.
            Used to annotate insertion gene assignments.

    Returns:
        pandas.DataFrame: Annotated insertions with extra metadata.

    """

    if isinstance(gtf, str):
        gtf = GtfFile(gtf)

    # Look-up genes in GTF file.
    genes = GtfFrame.from_records(gtf.fetch(filters={'feature': 'gene'}))
    genes.set_index('gene_id', drop=False, inplace=True)

    # Generate metadata.
    metadata = pd.DataFrame.from_records(
        (_annotate_insertion(ins, genes.ix[ins['gene_id']])
         for _, ins in insertions.iterrows()
         if ins['gene_id'] in genes.index))

    # Re-order columns.
    extra_cols = set(metadata.columns) - {'id', 'gene_id'}
    metadata = metadata[['id', 'gene_id'] + sorted(extra_cols)]

    return pd.merge(insertions, metadata, on=['id', 'gene_id'], how='left')


def _annotate_insertion(insertion, feature):
    """Annotates a given insertion/feature combination."""

    return {
        'id': insertion['id'],
        'gene_id': feature['gene_id'],
        'distance': feature_distance(insertion, feature),
        'orientation': feature_orientation(insertion, feature)
    }


def feature_distance(insertion, feature):
    """Calculates the genomic distance between an insertion and a feature.

    Args:
        insertion (pandas.Series): Insertion of interest. Assumed to have
            'position' and 'strand' values.
        feature (pandas.Series): Feature of interest. Assumed to have
            'start', 'end' and 'strand' values.

    Returns:
        int: Distance between insertion and feature.

    """

    feat_start, feat_end = feature['start'], feature['end']
    ins_location = insertion['position']

    if feat_start <= ins_location <= feat_end:
        dist = 0
    elif ins_location > feat_end:
        dist = ins_location - feat_end
    else:
        dist = ins_location - feat_start

    dist *= numeric_strand(feature['strand'])

    return dist


def feature_orientation(insertion, feature):
    """Determines the relative orientation of an insertion and a feature.

    Args:
        insertion (pandas.Series): Insertion of interest. Assumed to have
            a 'strand' value.
        feature (pandas.Series): Feature of interest. Assumed to have
            a 'strand' value.

    Returns:
        str: Returns 'sense' if features have the same orientation (i.e. are
            on the same strand), or 'antisense' if this is not the case.

    """

    ins_strand = numeric_strand(insertion['strand'])
    feat_strand = numeric_strand(feature['strand'])

    return 'sense' if ins_strand == feat_strand else 'antisense'
