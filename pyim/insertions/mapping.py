
import pandas


def map_insertions(hits, barcode_alignments):
    bc_lookup = barcode_alignments[['query_name', 'target_name']].set_index('query_name')
    bc_lookup.columns = ['barcode']

    hits_indexed = hits.set_index('query_name')
    put_insertions = hits_indexed.merge(bc_lookup, left_index=True, right_index=True, how='inner')

    return insertion_statistics(put_insertions)


def insertion_statistics(bc_hits):
    fwd_ins_stats = insertion_statistics_strand(bc_hits, '+', 'start', 'end')
    rev_ins_stats = insertion_statistics_strand(bc_hits, '-', 'end', 'start')
    return pandas.concat([fwd_ins_stats, rev_ins_stats], ignore_index=True)


def insertion_statistics_strand(bc_hits, strand, lp_field, grp_field):
    agg_funcs = { 'len': len, 'unique_len': lambda x: len(x.unique()) }

    str_ins = bc_hits.ix[bc_hits['strand'] == strand]
    str_ins_grp = str_ins.groupby(['barcode', 'contig_name', lp_field])
    str_ins_stats = str_ins_grp[grp_field].agg(agg_funcs)

    stats_names = ['barcode', 'chromosome', 'location', 'lp', 'unique_lp']

    str_ins_stats.reset_index(inplace=True)
    str_ins_stats.columns = stats_names
    str_ins_stats['strand'] = strand

    return str_ins_stats