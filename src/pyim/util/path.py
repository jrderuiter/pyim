from pathlib import Path


def build_path(file_path, suffix='', dir_=None, ext=None):
    file_path = Path(file_path)

    ext = ext or file_path.suffixes[-1]
    suffix = suffix + ext
    new_path = file_path.with_suffix(suffix)

    if dir_ is not None:
        new_path = Path(dir_) / new_path.name

    return new_path
