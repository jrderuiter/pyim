from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
from builtins import (ascii, bytes, chr, dict, filter, hex, input,
                      int, map, next, oct, open, pow, range, round,
                      str, super, zip)

import subprocess


def align(m1, index, output, m2=None, options=None, log=None):
    options = {} or options

    # Inject inputs into options.
    if m2 is None:
        options['-U'] = m1
    else:
        options['-1'] = m1
        options['-2'] = m2

    # Inject index and output.
    options['-x'] = index
    options['-S'] = output

    # Format into arguments.
    args = ['bowtie2'] + dict_to_args(options)

    if log is not None:
        with open(log, 'w') as log_file:
            subprocess.check_call(args, stderr=log_file)
    else:
        subprocess.check_call(args)

    return output


def dict_to_args(arg_dict):
    args = []

    for key, value in arg_dict.items():
        if type(value) == bool:
            if value:
                args.append(key)
        else:
            args.append(key)
            args.append(str(value))

    return args
