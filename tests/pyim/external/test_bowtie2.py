from pathlib import Path

import pytest

from pyim.external.bowtie2 import bowtie2, shell

# pylint: disable=redefined-outer-name


@pytest.fixture
def bowtie_args():
    """Basic arguments for bowtie2 function."""

    return {
        'in1_paths': [Path('/path/to/reads.R1.fastq')],
        'in2_paths': [Path('/path/to/reads.R2.fastq')],
        'output_path': Path('/path/to/output.bam'),
        'index_path': Path('/path/to/index'),
        'options': {'--threads': 10},
    }


def test_paired(mocker, bowtie_args):
    """Tests paired-end invocation of bowtie2."""

    mock = mocker.patch.object(shell, 'run_piped')
    bowtie2(**bowtie_args)

    expected_bt2 = ['bowtie2', '--threads', '10', '-1',
                    str(bowtie_args['in1_paths'][0]), '-2',
                    str(bowtie_args['in2_paths'][0]), '-x',
                    str(bowtie_args['index_path'])]
    expected_st = ['samtools', 'sort', '-o', str(bowtie_args['output_path']),
                   '-']

    mock.assert_called_with([expected_bt2, expected_st])


def test_single(mocker, bowtie_args):
    """Tests single-end invocation of bowtie2."""

    bowtie_args['in2_paths'] = None

    mock = mocker.patch.object(shell, 'run_piped')
    bowtie2(**bowtie_args)

    expected_bt2 = ['bowtie2', '--threads', '10', '-U',
                    str(bowtie_args['in1_paths'][0]), '-x',
                    str(bowtie_args['index_path'])]
    expected_st = ['samtools', 'sort', '-o', str(bowtie_args['output_path']),
                   '-']

    mock.assert_called_with([expected_bt2, expected_st])


def test_single_fa(mocker, bowtie_args):
    """Tests single-end invocation of bowtie2 with fasta file."""

    bowtie_args['in1_paths'] = [bowtie_args['in1_paths'][0].with_suffix('.fa')]
    bowtie_args['in2_paths'] = None

    mock = mocker.patch.object(shell, 'run_piped')
    bowtie2(**bowtie_args)

    expected_bt2 = ['bowtie2', '--threads', '10', '-U',
                    str(bowtie_args['in1_paths'][0]), '-f', '-x',
                    str(bowtie_args['index_path'])]
    expected_st = ['samtools', 'sort', '-o', str(bowtie_args['output_path']),
                   '-']

    mock.assert_called_with([expected_bt2, expected_st])
