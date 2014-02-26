
from pyim.io import read_fasta, write_fasta, Sequence
from pyim.util import setup_logging

from pyim.vector.sb import SBVectorAligner, SBAlignmentStats
from pyim.ggplot import ggplot_ipython, ggplot_ipython_svg
from pyim.genomic.blat import BlatAligner

from pyim.alignment.aligners import *
from pyim.alignment.filters import *

from IPython.core.display import display, display_svg


vector_file = '/Users/j.d.ruiter/Projects/NKI/IM/py-im-mapping/data/vec.fa'
barcode_file = '/Users/j.d.ruiter/Projects/NKI/IM/py-im-mapping/data/SB_barcodes.fa'

reference = '/Users/j.d.ruiter/References/Sanger_mm10/sequence/GRCm38_68.fa'

fasta_sanger1 = '/Users/j.d.ruiter/Work/20130826-sjors-im/readsets/1.TCA.454Reads.fna'
fasta_sanger2 = '/Users/j.d.ruiter/Work/20130826-sjors-im/readsets/2.TCA.454Reads.fna'
fasta_sb = '/Users/j.d.ruiter/Work/20130826-sjors-im/readsets/SBecad.TCA.454Reads.fna'
fasta_pool31 = '/Users/j.d.ruiter/Work/20130826-sjors-im/readsets/Pool3.1.TCA.454Reads.fna'

setup_logging()

vector_seqs = { s.name: s for s in read_fasta(vector_file) }
t7_seq = vector_seqs['T7']
sb_seq = vector_seqs['SB']


def show(x, width=400, height=300):
    display(ggplot_ipython(x, width=width, height=height))


def show_svg(x, width=7, height=5):
    display_svg(ggplot_ipython_svg(x, width=width, height=height))