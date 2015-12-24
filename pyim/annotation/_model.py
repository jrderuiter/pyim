
class Window(object):

    def __init__(self, start, end, reference=None, strand=None,
                 incl_left=True, incl_right=True, name=None):
        self.reference = reference
        self.start = start
        self.end = end
        self.strand = strand

        self.incl_left = incl_left
        self.incl_right = incl_right
        self.name = name

        if not incl_left or not incl_right:
            raise NotImplementedError()

    def apply(self, reference, location, strand=None):
        """Applies window to specific location and strand"""
        if strand is not None and self.strand is not None:
            strand = self.strand * strand
        else:
            strand = self.strand

        return Window(self.start + location, self.end + location,
                      reference, strand, self.incl_left,
                      self.incl_right, name=self.name)


# def apply_window(seqname, location, strand, window):
#     # TODO: Check strand logic!
#     start = location + (window.start * strand)
#     end = location + (window.end * strand)
#
#     if strand == -1:
#         start, end = end, start
#         incl_left, incl_right = window.incl_right, window.incl_left
#     else:
#         incl_left, incl_right = window.incl_left, window.incl_right
#
#     new_strand = strand * window.strand if window.strand is not None else None
#
#     return Window(seqname, start, end, new_strand, incl_left, incl_right)
