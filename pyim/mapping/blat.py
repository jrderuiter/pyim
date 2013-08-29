__author__ = 'j.d.ruiter'

import os
import subprocess
import pandas
from os import path

from pyim.io import read_psl, write_fasta


def blat_sequences(genomic_seqs, reference, output_dir='.', minscore=10, keep_fasta=True):
    fasta_file = path.join(output_dir, 'genomic.fna')
    write_fasta(genomic_seqs, fasta_file)

    psl_path = path.join(output_dir, 'genomic.psl')
    log_path = path.join(output_dir, 'blat.log')

    blat_cmd = "blat {reference} {fasta} -out=psl -minScore={minscore} {psl} 2> {log}"
    blat_cmd_frmt = blat_cmd.format(reference=reference, fasta=fasta_file,
                                    minscore=minscore, psl=psl_path, log=log_path)
    print blat_cmd_frmt

    subprocess.check_call(blat_cmd_frmt, shell=True)
    if not keep_fasta: os.unlink(fasta_file)

    return psl_path


def map_to_insertions(psl_file, bc_alignments):
    psl_frame = read_psl(psl_file)
    hits_frame = _filter_hits(psl_frame)

    bc_lookup = bc_alignments[['query_name', 'target_name']].set_index('query_name')
    bc_hits_frame = hits_frame.merge(bc_lookup, left_index=True, right_index=True, how='left')

    return insertion_statistics(bc_hits_frame)


def _filter_hits(psl_frame):    
    identity_mask = psl_frame['match'] >= psl_frame['q_size'] * 0.95

    masked = psl_frame[identity_mask]
    unique_reads = masked.drop(masked.index.get_duplicates())

    return unique_reads


def insertion_statistics(bc_hits):
    fwd_ins_stats = insertion_statistics_strand(bc_hits, '+', 't_start', 't_end')
    rev_ins_stats = insertion_statistics_strand(bc_hits, '-', 't_end', 't_start')
    return pandas.concat([fwd_ins_stats, rev_ins_stats], ignore_index=True)


def insertion_statistics_strand(bc_hits, strand, lp_field, grp_field):
    agg_funcs = { 'len': len, 'unique_len': lambda x: len(x.unique()) }

    str_ins = bc_hits.ix[bc_hits['strand'] == strand]
    str_ins_grp = str_ins.groupby(['target_name', 't_name', lp_field])
    str_ins_stats = str_ins_grp[grp_field].agg(agg_funcs)

    stats_names = ['barcode', 'chromosome', 'location', 'lp', 'unique_lp']

    str_ins_stats.reset_index(inplace=True)
    str_ins_stats.columns = stats_names
    str_ins_stats['strand'] = strand

    return str_ins_stats
