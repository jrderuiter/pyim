__author__ = 'j.d.ruiter'

import uuid, re, shutil, tempfile, pandas
import pandas.rpy.common as com
import rpy2.robjects as ro

from os import path
from rpy2.robjects.packages import importr
import rpy2.robjects.lib.ggplot2 as ggplot2


from IPython.core.display import Image
from IPython.display import SVG


grdevices = importr('grDevices')
r_print = ro.globalenv.get('print')

r_xlab = ro.globalenv.get('xlab')
r_ylab = ro.globalenv.get('ylab')
r_labs = ro.globalenv.get('labs')


def ggplot_ipython(plot, width=400, height=350, res=72):
    tmp_dir = tempfile.mkdtemp()
    fn = path.join(tmp_dir, '{uuid}.png').format(uuid = uuid.uuid4())

    grdevices.png(fn, width=width, height=height, res=res)
    plot.plot()
    grdevices.dev_off()

    im = Image(filename=fn)

    try: shutil.rmtree(tmp_dir)
    except OSError as exc:
        if exc.errno != 2: raise

    return im


def ggplot_ipython_svg(plot, width=2, height=1):
    tmp_dir = tempfile.mkdtemp()
    im_uuid = uuid.uuid4()

    fn = path.join(tmp_dir, '{uuid}.svg').format(uuid = im_uuid)

    grdevices.svg(fn, width=width, height=height)
    plot.plot()
    grdevices.dev_off()

    im = SVG(filename=fn)
    im._data = re.sub('#([a-zA-Z0-9-]+)', "#%s-\g<1>" % im_uuid.hex, im._data)
    im._data = re.sub('id="([a-zA-Z0-9-]+)"', 'id="%s-\g<1>"' % im_uuid.hex, im._data)

    try: shutil.rmtree(tmp_dir)
    except OSError as exc:
        if exc.errno != 2: raise

    return im


def ggplot_save(plot, filepath):
    grdevices.pdf(filepath)
    r_print(plot)
    grdevices.dev_off()


def _counts_to_rframe(counts):
    counts = pandas.DataFrame(counts, columns=['counts']).reset_index()
    return com.convert_to_r_dataframe(counts, strings_as_factors=True)


def ggplot_bar_hist(counts, title='', xlabel='', rotate=False, dense=False):
    r_counts = _counts_to_rframe(counts)

    pp = ggplot2.ggplot(r_counts) + \
         ggplot2.aes_string(x='index', y='counts') + \
         ggplot2.geom_bar(stat='identity', fill='#0072B2') + \
         ggplot2.scale_y_continuous('Counts')

    if dense: pp += ggplot2.geom_text(ggplot2.aes_string(label='counts'), angle=90, hjust=-0.1, size=2.5)
    else:     pp += ggplot2.geom_text(ggplot2.aes_string(label='counts'), vjust=-0.5, size=4)

    if title:  pp += ggplot2.ggtitle(title)
    if xlabel: pp += ggplot2.scale_x_discrete(xlabel)
    if rotate: pp += ggplot2.theme(**{'axis.text.x': ggplot2.element_text(angle=90, hjust=1)})

    return pp


def ggplot_stacked_bar(df, xvar, yvar, fillvar, xlabel='', ylabel='', title=''):
    r_frame = com.convert_to_r_dataframe(df, strings_as_factors=True)

    pp = ggplot2.ggplot(r_frame) + \
         ggplot2.aes_string(x='factor(%s)' % xvar, y=yvar, fill='factor(%s)' % fillvar) + \
         ggplot2.geom_bar(position='fill', stat='identity')

    if title: pp += ggplot2.ggtitle(title)
    if xlabel: pp += ggplot2.scale_x_discrete(xlabel)
    if ylabel: pp += ggplot2.scale_y_continuous(ylabel)

    return pp


def ggplot_hist(df, x='x', binwidth=None):
    if binwidth is None:
        binwidth = (df[x].max() - df[x].min()) / float(30)

    r_frame = com.convert_to_r_dataframe(df, strings_as_factors=True)
    pp = ggplot2.ggplot(r_frame) + \
         ggplot2.aes_string(x=x) + \
         ggplot2.geom_histogram(binwidth=binwidth)

    return pp


def ggplot_pie(counts):
    norm_counts = counts / float(counts.sum())
    r_counts = _counts_to_rframe(norm_counts)

    pp = ggplot2.ggplot(r_counts) + \
         ggplot2.aes_string(x='factor(1)', y='counts', fill='factor(index)') + \
         ggplot2.geom_bar(width=1, stat='identity') + ggplot2.coord_polar(theta='y') + \
         r_xlab('') + r_ylab('') + r_labs(fill='Response')

    return pp


