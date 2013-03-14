__author__ = 'j.d.ruiter'

import logging
import numpy as np
import pandas as pd


def map_to_barcodes(fasta_seqs, target_dict, aligner=None):
    logging.info('Mapping %d reads to barcodes with exact mapping', len(fasta_seqs))
    unmatched_seqs, exact_matched_map = barcode_exact_match(fasta_seqs, target_dict)

    if aligner is not None:
        logging.info('Mapping %d reads to barcodes with provided aligner', len(unmatched_seqs))
        unmatched_seqs, aligned_matched_map = barcode_alignment_match(unmatched_seqs, target_dict, aligner)
    else:
        aligned_matched_map = None

    return exact_matched_map, aligned_matched_map, unmatched_seqs


def barcode_exact_match(fasta_seqs, target_dict):


    seq_dict = {s.seqId: s for s in fasta_seqs}
    seq_ids = seq_dict.keys()
    seqs = [s.seq for s in seq_dict.values()]

    match_masks = dict()
    for tgt_id, tgt_seq in target_dict.items():
        match_masks[tgt_id] = [tgt_seq in s for s in seqs]

    match_frame = pd.DataFrame(match_masks, index=seq_ids)
    unmatched_list, matched_map = split_matches(match_frame)
    unmatched_seqs = [seq_dict[seq_id] for seq_id in unmatched_list]

    logging.info('- Exact matches: %d', len(matched_map))
    logging.info('- Unmatched: %d', len(unmatched_list))

    return unmatched_seqs, matched_map


def barcode_alignment_match(fasta_seqs, target_dict, aligner):
    alignments, scores = dict(), dict()
    for tgt_id, tgt_seq in target_dict.items():
        logging.info('- Aligning to %s', tgt_id)
        alignments[tgt_id] = aligner.align_target(fasta_seqs, tgt_seq)
        scores[tgt_id] = alignments[tgt_id]['score']
    score_frame = pd.DataFrame(scores)



    import pdb; pdb.set_trace()


def split_matches(match_frame):
    num_matches = match_frame.sum(axis=1)

    num_multiple = np.sum(num_matches > 1)
    if num_multiple > 0:
        logging.warning('Multiple exact matches were found for %d sequences, discarding', num_multiple)

    unmatched_list = list(match_frame[num_matches == 0].index)

    single_matches = match_frame[num_matches == 1]
    match_indices = single_matches.apply(lambda x: np.nonzero(x)[0][0], axis=1)

    target_ids = np.array(match_frame.columns)
    match_target_ids = target_ids[match_indices]
    match_series = pd.Series(match_target_ids, index=match_indices.index)

    return unmatched_list, match_series
