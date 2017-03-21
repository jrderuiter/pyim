from pathlib import Path

import pytest

from pyim.external.cutadapt import cutadapt, shell

# pylint: disable=redefined-outer-name


@pytest.fixture
def cutadapt_args():
    """Basic arguments for bowtie2 function."""

    return {
        'read_path': Path('/path/to/reads.R1.fastq'),
        'read2_path': Path('/path/to/reads.R2.fastq'),
        'out_path': Path('/path/to/output.R1.fastq'),
        'out2_path': Path('/path/to/output.R2.fastq'),
        'options': {'-m': 10},
    }


def test_paired(mocker, cutadapt_args):
    """Tests paired-end invocation of cutadapt."""

    mock = mocker.patch.object(shell, 'run')
    cutadapt(**cutadapt_args)

    expected = ['cutadapt', '-m', '10',
                '-o', str(cutadapt_args['out_path']),
                '-p', str(cutadapt_args['out2_path']),
                str(cutadapt_args['read_path']),
                str(cutadapt_args['read2_path'])] # yapf: disable
    mock.assert_called_with(expected)


def test_single(mocker, cutadapt_args):
    """Tests single-end invocation of cutadapt."""

    cutadapt_args['read2_path'] = None
    cutadapt_args['out2_path'] = None

    mock = mocker.patch.object(shell, 'run')
    cutadapt(**cutadapt_args)

    expected = ['cutadapt', '-m', '10',
                '-o', str(cutadapt_args['out_path']),
                str(cutadapt_args['read_path'])] # yapf: disable
    mock.assert_called_with(expected)
