__author__ = 'j.d.ruiter'
<<<<<<< HEAD


import pandas
from pyim.plot.hist import ggplot_hist


def hist_alignment_quality(full_alignment, extra_info):
    (vec_alignments, vec_unmapped), (bc_alignments, bc_unmapped) = extra_info

    t7_unmapped_ids = set([s.name for s in vec_unmapped['T7']])
    sb_unmapped_ids = set([s.name for s in vec_unmapped['SB']])
    bc_unmapped_ids = set([s.name for s in bc_unmapped])

    t7_sb_unmapped = t7_unmapped_ids & sb_unmapped_ids

    counts = pandas.Series({ 'Mapped': len(full_alignment),
                             'No SB': len(sb_unmapped_ids),
                             'No T7': len(t7_unmapped_ids),
                             'No BC': len(bc_unmapped_ids),
                             'No SB and T7': len(t7_sb_unmapped) })

    plot = ggplot_hist(counts, xlabel='Alignment status')

    return plot, counts
=======
>>>>>>> 0423cd6035e4f32a08eee4c112897caef8c3a66f
