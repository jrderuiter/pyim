

class Annotator(object):

    def annotate_by_gene(self, insertions):
        raise NotImplementedError

    def annotate_by_transcript(self, insertions):
        raise NotImplementedError
