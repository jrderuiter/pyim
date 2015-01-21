from collections import namedtuple


VectorAlignment = namedtuple(
    'VectorAlignment',
    ['q_name', 'q_start', 'q_end', 
     'v_name', 'v_start', 'v_end'])
