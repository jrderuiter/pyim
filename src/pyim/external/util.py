import subprocess


def run(arguments, stdout=None, stderr=None, check=True):
    stdout_ = _open_stdstream(stdout)
    stderr_ = _open_stdstream(stderr)

    try:
        process = subprocess.Popen(arguments, stdout=stdout_, stderr=stderr_)
        process.wait()
    finally:
        for std in [stdout_, stderr_]:
            _close_stdstream(std)

    # Check return code.
    if check and process.returncode != 0:
        raise ValueError('Process terminated with errorcode {}'
                         .format(process.returncode))

    return process


def run_piped(arguments_list, stdout=None, stderrs=None, check=True):
    if len(arguments_list) < 2:
        raise ValueError('At least two sets of arguments should be given')

    if stderrs is None:
        stderrs = [None] * len(arguments_list)

    # Handle processes 1 to n-1.
    processes = []
    stream_handles = []

    try:
        prev_out = None
        for arg_list, stderr in zip(arguments_list[:-1], stderrs[:-1]):
            # Setup processes.
            stderr_fh = _open_stdstream(stderr)
            stream_handles.append(stderr_fh)

            process = subprocess.Popen(
                arg_list,
                stdin=prev_out,
                stdout=subprocess.PIPE,
                stderr=stderr_fh)

            prev_out = process.stdout
            processes.append(process)

        # Handle final process.
        stdout_fh = _open_stdstream(stdout)
        stderr_fh = _open_stdstream(stderrs[-1])
        stream_handles += [stdout_fh, stderr_fh]

        process = subprocess.Popen(
            arguments_list[-1],
            stdout=stdout_fh,
            stderr=stderr_fh,
            stdin=prev_out)

        processes.append(process)

        # Allow pi to receive a SIGPIPE.
        for p in processes[:-1]:
            p.stdout.close()

        process.wait()

        # Check return codes.
        if check:
            if process.returncode != 0:
                raise ValueError('Process terminated with errorcode {}'
                                 .format(process.returncode))

    finally:
        # Close all file handles.
        for fh in stream_handles:
            _close_stdstream(fh)

    return processes


def _open_stdstream(file_path, mode='w'):
    if file_path is None:
        return subprocess.PIPE
    else:
        return file_path.open(mode)


def _close_stdstream(stdstream):
    if stdstream != subprocess.PIPE:
        stdstream.close()


def flatten_options(option_dict):
    """Flattens a dict of options into an argument list."""

    # Iterate over keys in lexical order, so that we have a
    # reproducible order of iteration (useful for tests).
    opt_names = sorted(option_dict.keys())

    # Flatten values.
    options = []
    for opt_name in opt_names:
        opt_value = option_dict[opt_name]

        if isinstance(opt_value, (tuple, list)):
            options += [opt_name] + [str(v) for v in opt_value]
        elif opt_value is True:
            options += [opt_name]
        elif not (opt_value is False or opt_value is None):
            options += [opt_name, str(opt_value)]

    return options
