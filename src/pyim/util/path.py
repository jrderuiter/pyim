import os

from pathlib import Path


def build_path(file_path, suffix='', dir_=None, ext=None):
    file_path = Path(file_path)

    try:
        ext = ext or file_path.suffixes[-1]
    except IndexError:
        ext = ''

    suffix = suffix + ext
    new_path = file_path.with_suffix(suffix)

    if dir_ is not None:
        new_path = Path(dir_) / new_path.name

    return new_path


def shorten_path(file_name, limit=40):
    """Shorten path for str to limit for logging."""

    name = os.path.split(str(file_name))[1]

    if len(name) > limit:
        return "%s~%s" % (name[:3], name[-(limit - 3):])
    else:
        return name
