from collections import namedtuple


GenomicAlignment = namedtuple(
    'GenomicAlignment',
    ['q_name', 'q_start', 'q_end', 'q_length',
     'r_name', 'r_start', 'r_end', 'r_strand'])
