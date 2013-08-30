__author__ = 'j.d.ruiter'

import uuid
import shutil
import tempfile
from os import path

import rpy2.robjects as ro
from rpy2.robjects.packages import importr

from IPython.core.display import Image


grdevices = importr('grDevices')
rprint = ro.globalenv.get('print')


def ggplot_ipython(plot, width=800, height=600, res=72):
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


def ggplot_save(plot, filepath):
    grdevices.pdf(filepath)
    rprint(plot)
    grdevices.dev_off()
