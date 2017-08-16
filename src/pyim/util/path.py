"""Utility functions for manipulating paths."""

import os

from contextlib import contextmanager
import tempfile
import shutil

from pathlib import Path


@contextmanager
def WorkDirectory(work_dir=None, keep=False):
    """Creates a work directory, which is optionally cleaned up."""

    if work_dir is None:
        work_dir, keep = Path(tempfile.mkdtemp()), False
    else:
        work_dir = Path(work_dir)
        work_dir.mkdir(exist_ok=True, parents=True)

    yield work_dir

    if not keep:
        shutil.rmtree(str(work_dir))


def shorten_path(file_name, limit=40):
    """Shorten path for str to limit for logging."""

    name = os.path.split(str(file_name))[1]

    if len(name) > limit:
        return "%s~%s" % (name[:3], name[-(limit - 3):])
    else:
        return name


def extract_suffix(file_path):
    """Extracts suffix from file path."""

    file_path = Path(file_path)

    if file_path.suffixes[-1] == '.gz':
        suffix = ''.join(file_path.suffixes[-2:])
    else:
        suffix = file_path.suffixes[-1]
    return suffix
