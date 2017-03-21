"""Utility module for running with external commands."""

import subprocess


def run(args, stdout=None, stderr=None, check=True):
    """Runs command for given arguments.

    Compared to the subprocess.run command, this function adds
    automatic opening/closing of any file paths passed to stdout/stderr.

    Parameters
    ----------
    args : List[str]
        Arguments for launching the process, passed as a list.
    stdout and stderr: Union[Path, int]
        These specify the executed programs' standard
        input, standard output and standard error file handles, respectively.
    check : bool
        Whether to check the returncode of the process.

    Returns
    -------
    subprocess.Popen
        Handle to completed process.

    """

    stdout_ = _open_stdstream(stdout)
    stderr_ = _open_stdstream(stderr)

    try:
        process = subprocess.Popen(args, stdout=stdout_, stderr=stderr_)
        process.wait()
    finally:
        for std in [stdout_, stderr_]:
            _close_stdstream(std)

    # Check return code.
    if check and process.returncode != 0:
        stderr_msg = process.stderr.read().decode()
        raise ValueError('Process terminated with errorcode {}\n\n'
                         'Output from stderr:\n\n'
                         .format(process.returncode) + stderr_msg)

    return process


def _open_stdstream(file_path, mode='w'):
    if file_path is None:
        return subprocess.PIPE
    else:
        return file_path.open(mode)


def _close_stdstream(stdstream):
    if stdstream != subprocess.PIPE:
        stdstream.close()


def run_piped(args_list, stdout=None, stderrs=None, check=True):
    """Runs piped commands for given argument lists."""

    if len(args_list) < 2:
        raise ValueError('At least two sets of arguments should be given')

    if stderrs is None:
        stderrs = [None] * len(args_list)

    # Handle processes 1 to n-1.
    processes = []
    stream_handles = []

    try:
        prev_out = None
        for arg_list, stderr in zip(args_list[:-1], stderrs[:-1]):
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

        final_process = subprocess.Popen(
            args_list[-1], stdout=stdout_fh, stderr=stderr_fh, stdin=prev_out)

        processes.append(final_process)

        # Allow pi to receive a SIGPIPE.
        for process in processes[:-1]:
            process.stdout.close()

        final_process.wait()

        # Check return codes.
        if check:
            if final_process.returncode != 0:
                stderr_msg = final_process.stderr.read().decode()
                raise ValueError('Process terminated with errorcode {}\n\n'
                                 'Output from stderr:\n\n'.format(
                                     final_process.returncode) + stderr_msg)

    finally:
        # Close all file handles.
        for handle in stream_handles:
            _close_stdstream(handle)

    return processes


def flatten_arguments(arg_dict):
    """Flattens a dict of options into an argument list.

    Parameters
    ----------
    arg_dict : Dict[Str, Any]
        Dictionary of arguments. Keys should be strings, values may be
        lists or tuples (for multiple values), booleans (for flags)
        or any other value that can be converted to a string.

    Returns
    -------
    List[str]
        List of flattened arguments.

    """

    # Iterate over keys in lexical order, so that we have a
    # reproducible order of iteration (useful for tests).
    opt_names = sorted(arg_dict.keys())

    # Flatten values.
    options = []
    for opt_name in opt_names:
        opt_value = arg_dict[opt_name]

        if isinstance(opt_value, (tuple, list)):
            options += [opt_name] + [str(v) for v in opt_value]
        elif opt_value is True:
            options += [opt_name]
        elif not (opt_value is False or opt_value is None):
            options += [opt_name, str(opt_value)]

    return options
