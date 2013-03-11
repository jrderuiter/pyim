
from os import path
import pandas as pd, os


class AlignmentCache(object):

    def __init__(self, fasta_path, cache_dir=None):
        if cache_dir is None:
            fasta_dir = path.dirname(fasta_path)
            cache_dir = path.join(fasta_dir, '_cache')

        fasta_name, _ = path.splitext(path.basename(fasta_path))

        self.fasta_path = fasta_path
        self.fasta_name = fasta_name
        self.cache_dir = cache_dir

    def add(self, target_name, alignment):
        file_path = self._cache_path(target_name)
        self._write(file_path, alignment)

    def get(self, target_name):
        file_path = self._cache_path(target_name)
        if path.exists(file_path):
            alignment = self._read(file_path)
        else:
            alignment = None
        return alignment

    def _cache_path(self, target_name):
        cache_name = '%s_%s.align.df' % (self.fasta_name, target_name)
        cache_path = path.join(self.cache_dir, cache_name)
        return cache_path

    def _write(self, file_path, alignment):
        target_dir = path.dirname(file_path)
        if not path.exists(target_dir):
            os.makedirs(target_dir)
        alignment.save(file_path)

    def _read(self, file_path):
        return pd.DataFrame.load(file_path)

    def __getitem__(self, item):
        return self.get(item)

    def __contains__(self, target_name):
        return path.exists(self._cache_path(target_name))
