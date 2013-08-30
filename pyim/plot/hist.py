__author__ = 'j.d.ruiter'

import pandas
import pandas.rpy.common as com

import rpy2.robjects.lib.ggplot2 as ggplot2


def ggplot_hist(counts, title='', xlabel='', rotate=False, dense=False):
    counts = pandas.DataFrame(counts, columns=['counts']).reset_index()
    r_counts = com.convert_to_r_dataframe(counts, strings_as_factors=True)

    pp = ggplot2.ggplot(r_counts) + \
         ggplot2.aes_string(x='index', y='counts') + \
         ggplot2.geom_bar(stat='identity', fill='#0072B2') + \
         ggplot2.scale_y_continuous('Counts')

    if dense:
        pp = pp + ggplot2.geom_text(ggplot2.aes_string(label='counts'), angle=90, hjust=-0.1, size=2.5)
    else:
        pp = pp + ggplot2.geom_text(ggplot2.aes_string(label='counts'), vjust=-0.5, size=4)

    if title: pp = pp + ggplot2.ggtitle(title)
    if xlabel: pp = pp + ggplot2.scale_x_discrete(xlabel)
    if rotate: pp = pp + ggplot2.theme(**{'axis.text.x': ggplot2.element_text(angle=90, hjust=1)})

    return pp


def hist_alignment_type(alignments, unmapped=None, title=''):
    counts = alignments['type'].value_counts()

    if unmapped is not None:
        counts = counts.append(pandas.Series(len(unmapped), index=['unmapped']))

    plot = ggplot_hist(counts, title, xlabel='Alignment type')

    return plot, counts


def hist_barcode_alignment(bc_alignments, bc_unmapped=None, title=''):
    bc_alignments = pandas.concat(bc_alignments.values(), ignore_index=True)

    bc_counts = bc_alignments['target_name'].value_counts()
    bc_counts = bc_counts.iloc[bc_counts.index.argsort()]

    if bc_unmapped is not None:
        bc_counts = bc_counts.append(pandas.Series(len(bc_unmapped), index=['unmapped']))

    plot = ggplot_hist(bc_counts, title, rotate=True, dense=True)

    return plot, bc_counts
