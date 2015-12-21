import collections


ExtractResult = collections.namedtuple(
    'ExtractResult', ['genomic_sequence', 'barcode', 'status'])


Insertion = collections.namedtuple(
    'Insertion', ['id', 'seqname', 'location',
                  'strand', 'sample', 'metadata'])
