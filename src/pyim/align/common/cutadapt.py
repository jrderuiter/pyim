import subprocess


def cutadapt(input_path, output_path, options):
    cmdline_args = _build_arguments(input_path, output_path, options)
    print(cmdline_args)
    #check_call(cmdline_args)


def cutadapt_piped(input_path, output_path, options_list):
    raise NotImplementedError()


def _build_arguments(input_path, output_path, options):
    """Builds argument list for cutadapt."""

    cmdline_opts = flatten_options(options)
    return (['cutadapt'] + cmdline_opts +
            ['-o', str(output_path), str(input_path)]) # yapf: disable


def _run(arguments, stdout=None, stderr=None, *args, **kwargs):
    stdout_ = _open_output(stdout)
    stderr_ = _open_output(stderr)

    try:
        process = subprocess.Popen(
            arguments, stdout=stdout, stderr=stderr, *args, **kwargs)
    finally:
        _close_output(stdout_)
        _close_output(stderr_)

    return process.returncode


def _run_piped(arguments_list, stdout=None, stderrs=None):
    if len(arguments_list) < 2:
        raise ValueError('At least two sets of arguments should be given')

    if stderrs is None:
        stderrs = [None] * len(arguments_list)

    # Handle processes 1-(n-1).
    processes = []
    file_handles = []

    try:
        prev_out = None
        for arg_list, stderr in list(zip(arguments_list, stderrs))[:-1]:
            # Setup processes.
            stderr_fh = _open_output(stderr)
            process = subprocess.Popen(
                arg_list,
                stdin=prev_out,
                stdout=subprocess.PIPE,
                stderr=stderr_fh)

            prev_out = process.stdout

            processes.append(process)
            file_handles.append(stderr_fh)

        # Handle final process.
        stdout_fh = _open_output(stdout)
        stderr_fh = _open_output(stderrs[-1])
        process = subprocess.Popen(
            arguments_list[-1],
            stdout=stdout_fh,
            stderr=stderr_fh,
            stdin=prev_out)

        processes.append(process)
        file_handles += [stderr_fh, stdout_fh]

        # Allow pi to receive a SIGPIPE.
        for p in processes[:-1]:
            p.stdout.close()

        process.wait()

    finally:
        # Close all file handles.
        for fh in file_handles:
            _close_output(fh)

    return process.returncode


def _open_output(file_path, mode='w'):
    if file_path is None:
        return None
    else:
        return file_path.open(mode)


def _close_output(file_path):
    if file_path is not None:
        file_path.close()


def flatten_options(option_dict):
    """Flattens a dict of options into an argument list."""

    options = []
    for opt_name, opt_values in option_dict.items():
        options += [opt_name] + list(opt_values)
    return options
