from collections import OrderedDict

from pyim.external import util


class TestRun(object):
    """Unit tests for the run function."""

    def test_simple(self):
        """Tests a simple command."""
        process = util.run(['echo', 'test'])

        assert process.returncode == 0
        assert process.stdout.read() == b'test\n'

    def test_stdout(self, tmpdir):
        """Tests a simple command with redirection to stdout."""

        log_path = tmpdir / 'test.log'
        process = util.run(['echo', 'test'], stdout=log_path)

        assert process.returncode == 0
        assert process.stdout is None

        with log_path.open('rb') as log_file:
            assert log_file.read() == b'test\n'


class TestRunPiped(object):
    """Unit tests for the run_piped function."""

    def test_simple(self):
        """Tests a simple piped command."""

        args = [['echo', 'test'], ['sed', 's/est/esting/g']]
        processes = util.run_piped(args)

        assert processes[-1].returncode == 0
        assert processes[-1].stdout.read() == b'testing\n'

    def test_stdout(self, tmpdir):
        """Tests a simple piped command with redirection to stdout."""

        log_path = tmpdir / 'test.log'
        args = [['echo', 'test'], ['sed', 's/est/esting/g']]
        processes = util.run_piped(args, stdout=log_path)

        assert processes[-1].returncode == 0
        assert processes[-1].stdout is None

        with log_path.open('rb') as log_file:
            assert log_file.read() == b'testing\n'


class TestFlattenArguments(object):
    """Unit tests for the flatten_arguments function."""

    def test_simple(self):
        """Tests a basic set of options."""

        options = OrderedDict([('--opt1', 'a'), ('--opt2', 'b')])
        result = util.flatten_arguments(options)

        assert result == ['--opt1', 'a', '--opt2', 'b']

    def test_non_str_type(self):
        """Tests options with a non-string type."""

        options = {'--opt1': 1}
        assert util.flatten_arguments(options) == ['--opt1', '1']

    def test_list(self):
        """Tests options with a list value."""

        options = OrderedDict([('--opt1', ['a', 'b']), ('--opt2', 'c')])
        result = util.flatten_arguments(options)
        assert result == ['--opt1', 'a', 'b', '--opt2', 'c']

    def test_flag_true(self):
        """Tests options with positive flag."""

        options = {'--opt1': True}
        assert util.flatten_arguments(options) == ['--opt1']

    def test_flag_false(self):
        """Tests options with negative flag."""

        options = {'--opt1': False}
        assert util.flatten_arguments(options) == []

    def test_flag_none(self):
        """Tests options with negative (None) flag."""

        options = {'--opt1': None}
        assert util.flatten_arguments(options) == []
