

class Annotator(object):

    def __init__(self, **kwargs):
        super(Annotator, self).__init__()

    def annotate_by_gene(self, insertions):
        raise NotImplementedError

    def annotate_by_transcript(self, insertions):
        raise NotImplementedError

    @classmethod
    def configure_argparser(cls, parser):
        return parser
