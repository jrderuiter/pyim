
import pandas


def hits_to_insertions(genomic_hits, barcode_mapping):
    barcoded_hits = genomic_hits
    barcoded_hits['barcode'] = genomic_hits['query_name'].map(barcode_mapping)

    # Collect insertions per strand and merge.
    insertions_fwd = _hits_to_insertions_strand(barcoded_hits, '+', 'start', 'end')
    insertions_rev = _hits_to_insertions_strand(barcoded_hits, '-', 'end', 'start')
    insertions = pandas.concat([insertions_fwd, insertions_rev], ignore_index=True)

    # Invert strand for insertion as the orientation of the sequence read
    # is the inverse orientation of the actual transposon insertion.
    insertions['strand'] = insertions['strand'].map({ '+': '-', '-': '+' })

    return insertions


def _hits_to_insertions_strand(barcoded_hits, strand, sb_field, t7_field):
    # Select hits on the specified strands. These hits are grouped by
    # barcode (to group them by sample) and by the position of the
    # SB side of the read (to group them by insertion in the sample).
    strand_hits = barcoded_hits.ix[barcoded_hits['strand'] == strand]
    strand_hits_grouped = strand_hits.groupby(['barcode', 'contig_name', sb_field])

    # Count the number of (unique) T7 ligation points
    agg_funcs = { 'len': len, 'unique_len': lambda x: len(x.unique()) }
    hit_group_stats = strand_hits_grouped[t7_field].agg(agg_funcs)

    # Create insertion frame from statistics
    insertions = hit_group_stats.reset_index()
    insertions.columns = ['barcode', 'chromosome', 'location', 'lp', 'unique_lp']
    insertions['strand'] = strand

    return insertions
