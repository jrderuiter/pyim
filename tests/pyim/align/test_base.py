import pytest

from pyim.align import base
from pyim.align.external.bowtie2 import shell as bt_shell
from pyim.align.external.cutadapt import shell as ct_shell
from pyim.model import InsertionSet


@pytest.fixture
def trim_kws():
    """Base kws for trim_transposon and trim_contaminants."""

    return {
        'input_paths': ('input/reads.R1.fastq.gz', 'input/reads.R2.fastq.gz'),
        'output_paths': ('out/trimmed.R1.fastq.gz', 'out/trimmed.R2.fastq.gz'),
        'sequences': {
            '-G': 'TGTATGTAAACTTCCGACTTCAACTG'
        },
        'threads': 4,
        'verbose': False
    }


class TestTrimTransposon(object):
    """Tests for the trim_transposon function."""

    @pytest.fixture
    def paired_kws(self, trim_kws):
        trim_kws['sequences'] = {'-G': 'TGTATGTAAACTTCCGACTTCAACTG'}
        return trim_kws

    @pytest.fixture
    def single_kws(self, paired_kws):
        paired_kws['input_paths'] = (paired_kws['input_paths'][0], )
        paired_kws['output_paths'] = (paired_kws['output_paths'][0], )
        return paired_kws

    def test_paired(self, paired_kws, mocker):
        """Test paired example."""

        mock_run = mocker.patch.object(ct_shell, 'run')
        base.trim_transposon(**paired_kws)

        expected_args = [
            'cutadapt', '--discard-untrimmed', '--pair-filter=both', '-G',
            paired_kws['sequences']['-G'], '-j',
            str(paired_kws['threads']), '-o', paired_kws['output_paths'][0],
            '-p', paired_kws['output_paths'][1], paired_kws['input_paths'][0],
            paired_kws['input_paths'][1]
        ]
        mock_run.assert_called_with(expected_args)

    def test_single(self, single_kws, mocker):
        """Test single-end example."""

        mock_run = mocker.patch.object(ct_shell, 'run')
        base.trim_transposon(**single_kws)

        expected_args = [
            'cutadapt', '--discard-untrimmed', '-G',
            single_kws['sequences']['-G'], '-j',
            str(single_kws['threads']), '-o', single_kws['output_paths'][0],
            single_kws['input_paths'][0]
        ]
        mock_run.assert_called_with(expected_args)

    def test_mismatched_inputs(self, paired_kws, mocker):
        """Test mismatched inputs."""

        paired_kws['output_paths'] = (paired_kws['output_paths'][0], )
        mocker.patch.object(ct_shell, 'run')

        with pytest.raises(ValueError):
            base.trim_transposon(**paired_kws)


class TestTrimContaminants(object):
    """Tests for the trim_contaminants function."""

    @pytest.fixture
    def paired_kws(self, trim_kws):
        trim_kws['sequences'] = {
            '-a': 'CTGTCTCTTATA',
            '-A': 'CTGTCTCTTATA',
        }
        return trim_kws

    @pytest.fixture
    def single_kws(self, paired_kws):
        paired_kws['input_paths'] = (paired_kws['input_paths'][0], )
        paired_kws['output_paths'] = (paired_kws['output_paths'][0], )
        return paired_kws

    def test_paired(self, paired_kws, mocker):
        """Test paired example."""

        mock_run = mocker.patch.object(ct_shell, 'run')
        base.trim_contaminants(**paired_kws)

        expected_args = [
            'cutadapt', '-A', paired_kws['sequences']['-A'], '-a',
            paired_kws['sequences']['-a'], '-j',
            str(paired_kws['threads']), '-o', paired_kws['output_paths'][0],
            '-p', paired_kws['output_paths'][1], paired_kws['input_paths'][0],
            paired_kws['input_paths'][1]
        ]
        mock_run.assert_called_with(expected_args)

    def test_single(self, single_kws, mocker):
        """Test single-end example."""

        mock_run = mocker.patch.object(ct_shell, 'run')
        base.trim_contaminants(**single_kws)

        expected_args = [
            'cutadapt', '-A', single_kws['sequences']['-A'], '-a',
            single_kws['sequences']['-a'], '-j',
            str(single_kws['threads']), '-o', single_kws['output_paths'][0],
            single_kws['input_paths'][0]
        ]
        mock_run.assert_called_with(expected_args)

    def test_mismatched_inputs(self, paired_kws, mocker):
        """Test mismatched inputs."""

        paired_kws['output_paths'] = (paired_kws['output_paths'][0], )
        mocker.patch.object(ct_shell, 'run')

        with pytest.raises(ValueError):
            base.trim_transposon(**paired_kws)


class TestAlignReads(object):
    """Tests for align_reads function."""

    @pytest.fixture
    def paired_kws(self):
        return {
            'read_paths': ('input/reads.R1.fastq.gz',
                           'input/reads.R2.fastq.gz'),
            'index_path': '/path/to/index',
            'output_path': 'out/aligned.bam',
            'threads': 4,
            'verbose': False
        }

    @pytest.fixture
    def single_kws(self, paired_kws):
        paired_kws['read_paths'] = (paired_kws['read_paths'][0], )
        return paired_kws

    def test_paired(self, paired_kws, mocker):
        """Test paired example."""

        mock_run = mocker.patch.object(bt_shell, 'run_piped')
        base.align_reads(**paired_kws)

        bowtie_expected = [
            'bowtie2', '--threads', '4', '-1', 'input/reads.R1.fastq.gz', '-2',
            'input/reads.R2.fastq.gz', '-x', '/path/to/index'
        ]
        sort_expected = ['samtools', 'sort', '-o', 'out/aligned.bam', '-']

        mock_run.assert_called_with([bowtie_expected, sort_expected])

    def test_single(self, single_kws, mocker):
        """Test single example."""

        mock_run = mocker.patch.object(bt_shell, 'run_piped')
        base.align_reads(**single_kws)

        bowtie_expected = [
            'bowtie2', '--threads', '4', '-U', 'input/reads.R1.fastq.gz', '-x',
            '/path/to/index'
        ]
        sort_expected = ['samtools', 'sort', '-o', 'out/aligned.bam', '-']

        mock_run.assert_called_with([bowtie_expected, sort_expected])


class TestExtractInsertions(object):
    """Tests for the extract_insertions function."""

    @staticmethod
    def _position_for_mates(mate1, mate2):
        """Returns transposon/linker positions for given mates."""

        ref = mate1.reference_name

        if mate1.is_reverse:
            transposon_pos = mate2.reference_start
            linker_pos = mate1.reference_end
            strand = 1
        else:
            transposon_pos = mate2.reference_end
            linker_pos = mate1.reference_start
            strand = -1

        return (ref, transposon_pos, strand), linker_pos

    def test_paired(self):
        """Tests extraction of insertions from paired-end alignment."""
        bam_path = pytest.helpers.data_path('tagmap/subset.bam')

        insertions = InsertionSet.from_tuples(
            base.extract_insertions(
                bam_path,
                position_func=self._position_for_mates,
                sample_func=lambda m1, m2: 'sample',
                min_mapq=23,
                merge_dist=10,
                min_support=2,
                paired=True))

        assert len(insertions) == 7

        # Check IDs are numbered consecutively.
        assert set(insertions['id'] == {
            'sample.INS_{}'.format(i + 1)
            for i in range(7)
        })

        # Sort by positions and check values.
        insertions = insertions.sort_values(by=['chromosome', 'position'])

        assert all(
            insertions['chromosome'] == ['10', '10', '15', '2', '4', '7', '8'])

        assert all(insertions['position'] == [
            70280823, 108240765, 50807519, 59687915, 133730139, 130262794,
            124583415
        ])

        assert all(insertions['support'] == [16, 4, 4, 4, 6, 3, 2])
        assert all(insertions['strand'] == [-1, -1, 1, -1, 1, -1, -1])
