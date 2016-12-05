from pathlib import Path

import pytest

from pyim.util import path

# pylint: disable=redefined-outer-name


@pytest.fixture
def test_path():
    return Path('path/to/file.txt')


class TestBuildPath(object):
    def test_suffix(self, test_path):
        """Tests addition of suffix to path."""

        result = path.build_path(test_path, suffix='.filt')
        assert str(result) == 'path/to/file.filt.txt'

    def test_suffix_missing_ext(self, test_path):
        """Tests addition of suffix to path that is missing an ext."""

        result = path.build_path(test_path.with_suffix(''), suffix='.filt')
        assert str(result) == 'path/to/file.filt'

    def test_ext(self, test_path):
        """Tests replacing ext of path."""

        result = path.build_path(test_path, ext='.log')
        assert str(result) == 'path/to/file.log'

    def test_ext_missing_ext(self, test_path):
        """Tests replacing ext of path that is missing an ext."""

        result = path.build_path(test_path.with_suffix(''), ext='.log')
        assert str(result) == 'path/to/file.log'

    def test_dir(self, test_path):
        """Tests changing directory of path."""

        result = path.build_path(test_path, dir_='other')
        assert str(result) == 'other/file.txt'

    def test_nop(self, test_path):
        """Tests calling of build_path with no real arguments."""

        assert str(test_path) == str(path.build_path(test_path))

    def test_str_path(self, test_path):
        """Tests calling of build_path with str path."""

        assert str(test_path) == str(path.build_path(str(test_path)))
