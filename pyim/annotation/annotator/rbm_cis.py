from __future__ import absolute_import, division, print_function

#pylint: disable=wildcard-import,unused-wildcard-import,redefined-builtin
from builtins import *
#pylint: enable=wildcard-import,unused-wildcard-import,redefined-builtin

import logging
from os import path

import numpy as np
import pandas as pd

#pylint: disable=import-error
from .rbm import annotate_rbm, WINDOW_PRESETS as RBM_WINDOW_PRESETS
from ..metadata import add_metadata
from ..filtering import filter_blacklist, select_closest
#pylint: enable=import-error


def annotate_rbm_cis(insertions, cis_sites, gtf, window_preset=None,
                     window_sizes=None, collapse=False):
    """Assigns insertions to genes using the RBM approach via called CIS sites.

    Args:
        insertions (pandas.DataFrame): Insertions to annotate in DataFrame
            format. The frame is expected to contain at least the
            following columns: id, position, strand.
        cis_sites(pandas.DataFrame): Dataframe containing the CIS sites
            for the given insertions.
        gtf (str or GtfFile): Path to gtf file containing gene features.
            Alternatively, a GtfFile object may also be given instead of a path.
        window_preset (str): Preset to use for the RBM window sizes.
            Alternatively custom window sizes can be given using the
            *window_sizes* argument. Note that either *window_preset* or
            *window_sizes* must be provided.
        window_sizes (tuple[int]): Tuple of window sizes to use in the
            RBM mapping. Should specify four window sizes, for the following
            categories of insertions: upstream-sense, upstream-antisense,
            downstream-sense, downstream-antisense.

    Returns:
        tuple[pandas.DataFrame]: Returns two dataframes, the first
            containing the annotated insertion sites, the second containing
            the annotated CIS sites, which were used to annotate the insertions.
            Annotations are added as columns 'gene_id' and 'gene_name', which
            respectively contain the id and name of the annotated gene. An
            extra column 'window' indicates which of the RBM windows was
            used for the annotation.

    """

    if 'strand' not in cis_sites:
        # Add strand to cis sites if not present.
        cis_sites = _determine_cis_strand(cis_sites, insertions)

    # Annotate cis sites.
    cis_sites = cis_sites.rename(columns={'cis_id': 'id'})
    annotated_sites = annotate_rbm(cis_sites, gtf,
                                   window_preset=window_preset,
                                   window_sizes=window_sizes)

    # Extract and merge annotation with insertions.
    annotation = annotated_sites[['id', 'gene_id', 'gene_name']]
    annotation = annotation.rename(columns={'id': 'cis_id'})
    annotated_ins = pd.merge(insertions, annotation, on='cis_id', how='left')

    if collapse:
        # Collapse multiple insertion entries resulting from CIS annotation.
        annotated_ins.drop(['cis_id'], axis=1, inplace=True)
        annotated_ins.drop_duplicates(inplace=True)

    return annotated_ins, annotated_sites


def _determine_cis_strand(cis, cis_insertions, min_homogeneity=0.5):
    """Determines the strand for CIS sites with homogeneous insertions."""

    # Extract and clip strands at zero.
    ins_strands = cis_insertions[['cis_id', 'strand']].copy()
    ins_strands['strand'] = ins_strands['strand'].map({1: 1, -1: 0})

    # Calculate fwd/rev ratio for each cis.
    ratio = ins_strands.groupby('cis_id')['strand'].mean()

    # Determine closest strand and homogeneity.
    cis_strands = pd.DataFrame(
        {'strand': ratio.round().astype(int).map({1: 1, 0: -1}),
         'strand_homogeneity': np.maximum((1 - ratio), ratio)},
        columns=['strand', 'strand_homogeneity'])

    # Don't assign strand if low homogeneity.
    homogeneity_mask = cis_strands['strand_homogeneity'] < min_homogeneity
    cis_strands.ix[homogeneity_mask, 'strand'] = None

    return pd.merge(cis, cis_strands.reset_index())


def register(subparsers, name='rbm-cis'):
    """Registers the RBM-CIS annotator as a subparser."""

    parser = subparsers.add_parser(name, help=name + ' help')

    # Required arguments.
    parser.add_argument('input')
    parser.add_argument('output')
    parser.add_argument('--gtf', required=True)
    parser.add_argument('--cis_sites', required=True)

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--preset', choices=RBM_WINDOW_PRESETS.keys())
    group.add_argument('--window_sizes', nargs=4, type=int)

    # Optional arguments.
    parser.add_argument('--closest', default=False, action='store_true')
    parser.add_argument('--collapse', default=False, action='store_true')
    parser.add_argument('--blacklist', default=None, nargs='+')

    # Set main for dispatch.
    parser.set_defaults(main=main)

    return parser


def main(args):
    """Main function for the RBM-CIS annotator command-line tool."""

    # Read insertions and cis sites.
    insertions = pd.read_csv(args.input, sep='\t', dtype={'chrom': str})
    cis_sites = pd.read_csv(args.cis_sites, sep='\t', dtype={'chrom': str})

    logging.info('Read %d insertions and %d cis sites',
                 len(insertions), len(cis_sites))

    # Annotate insertions.
    logging.info('Annotating insertions')

    annotated_ins, annotated_sites = annotate_rbm_cis(
        insertions, cis_sites, args.gtf, window_preset=args.preset,
        window_sizes=args.window_sizes, collapse=args.collapse)

    # Add metadata to annotated insertions.
    logging.info('Adding annotation metadata')
    annotated_ins = add_metadata(annotated_ins, args.gtf)

    if args.blacklist is not None:
        logging.info('Filtering blacklisted genes')
        annotated_ins = filter_blacklist(annotated_ins, args.blacklist)
        annotated_sites = filter_blacklist(annotated_sites, args.blacklist)

    if args.closest:
        logging.info('Selecting closest insertions')
        annotated_ins = select_closest(annotated_ins)
        annotated_sites = add_metadata(annotated_sites, args.gtf)
        annotated_sites = select_closest(annotated_sites)

    # Write outputs.
    annotated_ins.to_csv(args.output, sep='\t', index=False)
    annotated_sites.to_csv(path.splitext(args.output)[0] + '.sites.txt',
                           sep='\t', index=False)
