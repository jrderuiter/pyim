__author__ = 'Julian'

import pandas

from pyim_common.util import reorder_columns


def alignment_frame(sam_file, alignments):
    dicts_ = (_aln_to_dict(sam_file, a) for a in alignments)
    frame = pandas.DataFrame(dicts_)
    frame = reorder_columns(frame, order=['name', 'seqname', 'start', 'end', 'strand'])
    return frame


def _aln_to_dict(sam_file, aln):
    return {
        'name': aln.query_name,
        'seqname': sam_file.getrname(aln.reference_id),
        'start': aln.pos,
        'end': aln.aend,
        'strand': -1 if aln.is_reverse else 1,
        'mapq': aln.mapping_quality
    }


def chunk_alignments(alignments, max_dist=5000):
    aln = next(alignments)

    # Setup initial chunk state.
    chunk = [aln]
    prev_seq = aln.reference_id
    prev_pos = aln.pos

    for aln in alignments:
        if aln.reference_id != prev_seq or (aln.pos - prev_pos) > max_dist:
            yield chunk

            chunk = [aln]
            prev_seq = aln.reference_id
        else:
            chunk.append(aln)

        prev_pos = aln.pos

    ## Don't forget our final yield.
    yield chunk
