from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
from builtins import (ascii, bytes, chr, dict, filter, hex, input,
                      int, map, next, oct, open, pow, range, round,
                      str, super, zip)

import os
import subprocess
from os import path


def align(m1, index, output, m2=None, options=None,
          log=None, bam_output=False):
    options = {} or options

    # Inject inputs into options.
    if m2 is None:
        options['-U'] = m1
    else:
        options['-1'] = m1
        options['-2'] = m2

    # Inject index and output.
    if not output.endswith('.sam'):
        if output.endswith('.bam'):
            output = output.replace('.bam', '.sam')
        else:
            output = output + '.sam'

    options['-x'] = index
    options['-S'] = output

    # Format into arguments.
    args = ['bowtie2'] + dict_to_args(options)

    if log is not None:
        with open(log, 'w') as log_file:
            subprocess.check_call(args, stderr=log_file)
    else:
        subprocess.check_call(args)

    # Convert to bam if needed.
    if bam_output:
        output = sam_to_bam(output, sort=True,
                            index=True, delete_sam=True)

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


def sam_to_bam(sam_path, bam_path=None, sort=False,
               index=False, delete_sam=False):
    if bam_path is None:
        # Default output name replaces .sam with .bam.
        bam_path = path.splitext(sam_path)[0] + '.bam'

    if sort:
        # Pipe bam into samtools sort for sorting.
        p1 = subprocess.Popen(['samtools', 'view', '-b', sam_path],
                              stdout=subprocess.PIPE)
        p2 = subprocess.Popen(['samtools', 'sort', '-o', bam_path, '-'],
                              stdin=p1.stdout)
        p1.stdout.close()
        p2.communicate()

        if index:
            # Indexing bam file if needed.
            subprocess.check_call(['samtools', 'index', bam_path])
    else:
        # Only convert to bam.
        subprocess.check_call(['samtools', 'view', '-b',
                               '-o', bam_path, sam_path])

    if delete_sam:
        # Delete original sam if requested.
        os.unlink(sam_path)

    return bam_path
