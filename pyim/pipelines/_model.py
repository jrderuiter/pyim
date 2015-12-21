import collections


ExtractResult = collections.namedtuple(
    'ExtractResult', ['genomic_sequence', 'barcode', 'status'])
