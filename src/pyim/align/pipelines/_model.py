
class ExtractResult(object):

    __slots__ = ('genomic_sequence', 'barcode', 'status')

    def __init__(self, genomic_sequence, barcode, status):
        super().__init__()
        self.genomic_sequence = genomic_sequence
        self.barcode = barcode
        self.status = status
