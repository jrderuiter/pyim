#!/usr/bin/env python

from __future__ import print_function, unicode_literals, \
                       absolute_import, division
import argparse
import logging
import pandas as pd
import numpy as np

from pyim.io.fasta import read_fasta
from pyim.io.psl import read_psl
from pyim.alignment.aligners.exact import ExactReadAligner
from pyim.alignment.aligners.inexact import WatermanAligner
from pyim.alignment.aligners.combined import CombinedReadAligner

logging.basicConfig(format='[%(asctime)s] %(levelname)-8s %(message)s', level=logging.INFO, datefmt='%H:%M:%S')

fasta_file = '/Volumes/Datastore/Julian/Datasets/sjors-im/integrations/Pool3/Pool3.1.TCA.454Reads.fna'
vector_file = './data/vec.fa'
barcode_file = './data/SBbarcodes.fa'
num_cores = 2
read_limit = 100
read_alignment = '/Volumes/Datastore/Julian/Datasets/sjors-im/integrations/Pool3/POOL3.1.psl'
min_bc_size = 1

def main(fasta_file, vector_file, barcode_file, num_cores, read_limit, read_alignment, min_bc_size):
    logging.info('Loading fasta sequences')
    fasta_seqs = _load_fasta(fasta_file, read_limit)

    # Prepare combined aligner
    aligner = CombinedReadAligner(ExactReadAligner(), WatermanAligner(num_processes=num_cores))

    # Align vector sequences and barcodes
    #import pdb; pdb.set_trace()
    vec_alignments, _ = _load_and_align_vectors(fasta_seqs, vector_file, aligner)
    sb_alignments = vec_alignments[vec_alignments['target_name'] == 'SB']

    bc_alignments, _ = _load_and_align_barcodes(fasta_seqs, barcode_file, aligner)
    valid_bc_alns, bc_groups = _filter_bc(bc_alignments, min_bc_size)

    # Load hits
    logging.info('Loading BLAT alignment')
    hits = read_psl(read_alignment, onlyPutative=False)

    joined = hits.add_prefix('hit_').join(sb_alignments.add_prefix('sb_'), how='left')
    joined = joined.join(bc_alignments.add_prefix('bc_'), how='left')

    # Apply filters
    sb_mask = _sb_mask(joined, prefix='sb_')
    hit_mask = _hit_mask(joined, prefix='hit_')
    dist_mask = _hit_sb_mask(joined, sb_prefix='sb_', hit_prefix='hit_')
    full_mask = sb_mask & hit_mask & dist_mask

    filtered_hits = joined.ix[full_mask]

    stats = []
    for name, group in filtered_hits.groupby('bc_target_name'):
        bc_hits_all = group
        bc_hits_putative = select_strongest(bc_hits_all, prefix='hit_')

        bc_strand = bc_hits_putative['hit_strand']
        bc_hits_forward = bc_hits_putative[bc_strand == '+']
        bc_hits_reverse = bc_hits_putative[bc_strand == '-']

        bc_forward_sub = bc_hits_forward[['hit_t_name', 'hit_t_start', 'hit_strand']]
        bc_reverse_sub = bc_hits_reverse[['hit_t_name', 'hit_t_end', 'hit_strand']]
        bc_reverse_sub.columns = ['name', 'start', 'strand']
        bc_forward_sub.columns = ['name', 'start', 'strand']
        bc_joined = pd.concat([bc_forward_sub, bc_reverse_sub])
        bc_joined_nodup = bc_joined.drop_duplicates()

        stats.append(dict(name=name, num_reads=len(group), num_mapping=len(group.index.unique()),
                          num_inserts=len(bc_hits_putative), num_uni_inserts=len(bc_joined_nodup)))
    stats = pd.DataFrame(stats)
    stats.reindex_axis(['name', 'num_reads', 'num_mapping', 'num_inserts', 'num_uni_inserts'], axis=1)



def _load_fasta(fasta_file, read_limit):
    logging.info('Loading reads from fasta input file')
    fastaSeqs = list(read_fasta(fasta_file))

    if read_limit:
        logging.info('- Limiting to %d reads' % read_limit)
        fastaSeqs = fastaSeqs[0:read_limit]

    return fastaSeqs

def _load_and_align_vectors(fasta_seqs, vector_file, aligner):
    vec_seq_dict = {s.seqId: s for s in read_fasta(vector_file)}
    vec_seqs = [vec_seq_dict[seq_id] for seq_id in ['T7', 'SB']]
    vec_alignments = aligner.align_targets(fasta_seqs, vec_seqs)
    return vec_alignments, vec_seqs

def _load_and_align_barcodes(fasta_seqs, barcode_file, aligner):
    barcode_seqs = list(read_fasta(barcode_file))[0:2]
    bc_alignments = aligner.align_targets(fasta_seqs, barcode_seqs, parallel_in_targets=True)
    return bc_alignments, barcode_seqs

def _filter_bc(bc_alignments, min_bc_size):
    bc_groups = bc_alignments.groupby('target_name')
    valid_bc_groups = filter(lambda x: len(x[1]) > min_bc_size, bc_groups)
    valid_frame = pd.concat([frame for _, frame in valid_bc_groups])
    return valid_frame, bc_groups

def _sb_mask(sb_frame, prefix=''):
    starts_at_begin = sb_frame[prefix + 'target_start'] == 0
    match_long_enough = sb_frame[prefix + 'target_end'] >= 23
    alignment_sufficient = sb_frame[prefix + 'score'] >= 75

    valid_mask = starts_at_begin & match_long_enough & \
                 alignment_sufficient
    return valid_mask

def select_strongest(hits, prefix=''):
    def selmax(group):
        match_scores = group[prefix + 'match']
        max_val = np.max(match_scores)
        return group.ix[match_scores == max_val]

    grps = hits.groupby(level=0)
    return grps.apply(selmax)

def _hit_mask(hit_frame, prefix=''):
    total = hit_frame[prefix + 'mismatch'] + hit_frame[prefix + 'match']
    hit_identity = hit_frame[prefix + 'match'] / total
    identity_mask = hit_identity >= 0.95
    return identity_mask

def _hit_sb_mask(joined_frame, sb_prefix='sb_', hit_prefix='hit_'):
    hit_dist = np.abs(joined_frame[hit_prefix + 'q_start'] - joined_frame[sb_prefix + 'query_end'])
    dist_mask = hit_dist <= 5
    return dist_mask

def _parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--reads', required=True)
    parser.add_argument('--vectors', default='vec.fa')
    parser.add_argument('--barcodes', default='data/SBbarcodes.fa')
    parser.add_argument('--read_alignment', required=True)
    parser.add_argument('--num_cores', type=int, default=1)
    parser.add_argument('--read_limit', type=int, default=None)

    return parser.parse_args()


if __name__ == '__main__':
    args = _parse_args()
    main(args.reads, args.vectors, args.barcodes, args.num_cores, args.read_limit, args.read_alignment, 1)