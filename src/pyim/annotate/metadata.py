from pathlib import Path

import toolz

from pyim.util.tabix import GtfFile, GtfFrame
from .util import numeric_strand


def add_metadata(insertions, reference_gtf):
    """Adds metadata to annotated insertions.

    Adds extra metadata to already annotated insertions. This metadata
    currently includes the following information: distance to the gene
    ('distance' column) and relative orientation ('orientation' column).

    Args:
        insertions (iterable[Insertion]): Annotated insertions for which metadata
            should be added. The frame is expected to contain at least the
            following columns: id, position, strand, gene_id.
        gtf (str, Path, GtfFile of GtfFrame): Gtf containing gene features.

    Returns:
        iterable[Insertion]: Annotated insertions with extra metadata.

    """

    if isinstance(reference_gtf, (str, Path)):
        reference_gtf = GtfFile(reference_gtf)

    if isinstance(reference_gtf, GtfFile):
        reference_gtf = GtfFrame.from_records(
            reference_gtf.fetch(filters={'feature': 'gene'}))

    # Look-up genes in GTF frame.
    genes = reference_gtf.get_region(filters={'feature': 'gene'})
    genes.set_index('gene_id', drop=False, inplace=True)

    for insertion in insertions:
        if 'gene_id' in insertion.metadata:

            # Add metadata for gene.
            gene = genes.ix[insertion.metadata['gene_id']]

            gene_metadata = {
                'gene_distance': feature_distance(insertion, gene),
                'gene_orientation': feature_orientation(insertion, gene)
            }

            new_metadata = toolz.merge(insertion.metadata, gene_metadata)
            yield insertion._replace(metadata=new_metadata)
        else:
            # Return original insertion.
            yield insertion


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
    ins_location = insertion.position

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

    ins_strand = numeric_strand(insertion.strand)
    feat_strand = numeric_strand(feature['strand'])

    return 'sense' if ins_strand == feat_strand else 'antisense'
