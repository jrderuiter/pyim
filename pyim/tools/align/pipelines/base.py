__author__ = 'Julian'

from pyim_common.util import flatten_list

from pyim.tools.align.util import chunk_alignments


class Pipeline(object):

    REQUIRED_OPTIONS = []

    def __init__(self, **kwargs):
        super(Pipeline, self).__init__()

    @classmethod
    def configure_argparser(cls, parser):
        return parser

    @staticmethod
    def map_insertions_chunked(sam_file, alignments, map_func):
        chunks = chunk_alignments(alignments)
        chunk_ins = [map_func(sam_file, c) for c in chunks]
        return flatten_list(chunk_ins)
