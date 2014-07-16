import logging

import pandas

from pyim.insertions.clustering import cluster_insertions
from pyim.insertions.mapping import map_insertions_to_samples

logger = logging.getLogger('PyIM.SB.Insertions')


def identify_insertions(alignments, barcode_sample_file, min_ulp):
    logger.info('Generating insertions from alignments')
    insertions = alignments_to_insertions(alignments)
    insertions = cluster_insertions(insertions, max_dist=10)

    logger.info('Restricting insertions to main chromosomes')
    main_chromosomes = list(map(str, range(1, 19+1))) + ['X', 'Y']
    insertions = insertions.ix[insertions['chromosome'].isin(main_chromosomes)]

    logger.info('Removing insertions with ulp < %d' % min_ulp)
    insertions = insertions.ix[insertions['unique_lp'] >= min_ulp]

    logger.info('Mapping insertions to samples')
    insertions['id'] = ['INS_%d' % ind for ind in range(1, len(insertions)+1)]
    insertions = map_insertions_to_samples(insertions, barcode_sample_file, drop_barcode=True)

    return insertions[['id', 'chromosome', 'location', 'strand', 'sample', 'lp', 'unique_lp']]


def alignments_to_insertions(alignments):
    insertions_fwd = _alignments_to_insertions_strand(alignments, '+', 'start', 'end')
    insertions_rev = _alignments_to_insertions_strand(alignments, '-', 'end', 'start')
    return pandas.concat([insertions_fwd, insertions_rev], ignore_index=True)


def _alignments_to_insertions_strand(alignments, strand, sb_field, t7_field):
    # Select hits on the specified strands. These hits are grouped by
    # barcode (to group them by sample) and by the position of the
    # SB side of the read (to group them by insertion in the sample).
    strand_hits = alignments.ix[alignments['strand'] == strand]
    strand_hits_grouped = strand_hits.groupby(['bc_name', 'contig_name', sb_field])

    # Count the number of (unique) T7 ligation points
    agg_funcs = {'len': len, 'unique_len': lambda x: len(x.unique())}
    hit_group_stats = strand_hits_grouped[t7_field].agg(agg_funcs)

    # Create insertion frame from statistics
    insertions = hit_group_stats.reset_index()
    insertions = insertions[['bc_name', 'contig_name', sb_field, 'len', 'unique_len']]
    insertions.columns = ['barcode', 'chromosome', 'location', 'lp', 'unique_lp']
    insertions['strand'] = strand

    return insertions
