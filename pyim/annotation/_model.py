
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

    def apply(self, reference, location, strand):
        """Applies window to specific location and strand"""

        # Determine start/end position.
        if strand == 1:
            start = location + self.start
            end = location + self.end

            incl_left = self.incl_left
            incl_right = self.incl_right
        elif strand == -1:
            start = location - self.end
            end = location - self.start

            incl_right = self.incl_left
            incl_left = self.incl_right
        else:
            raise ValueError('Unknown value for strand ({})'
                             .format(strand))

        # Determine new strand.
        if self.strand is not None:
            new_strand = self.strand * strand
        else:
            new_strand = None

        return Window(start, end, reference, new_strand,
                      incl_left, incl_right, name=self.name)
