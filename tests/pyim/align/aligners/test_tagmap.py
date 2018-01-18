from mock import ANY
from pathlib import Path
import shutil

import pytest

from pyim.align.aligners import TagmapAligner
from pyim.align.aligners import tagmap


class TestTagmapAligner(object):
    """Tests for the TagmapAligner class."""

    @pytest.fixture
    def subset_bam(self, tmpdir):
        """Copies example bam into expected location in tmpdir."""
        src_path = pytest.helpers.data_path('tagmap/subset.bam')
        dest_path = tmpdir / 'test.bam'
        shutil.copy(str(src_path), str(dest_path))
        return dest_path

    def test_example(self, subset_bam, mocker):
        """Test using example data."""

        # Mock 'external' functions called by aligner.
        mock_trim_tr = mocker.patch.object(tagmap, 'trim_transposon')
        mock_trim_nt = mocker.patch.object(tagmap, 'trim_contaminants')
        mock_align = mocker.patch.object(tagmap, 'align_reads')
        mock_extract = mocker.spy(tagmap, 'extract_insertions')

        # Build aligner.
        prefix = str(subset_bam).replace('.bam', '')

        aligner = TagmapAligner(
            transposon_seq='TGTATGTAAACTTCCGACTTCAACTG',
            bowtie_index_path=Path('/path/to/index'),
            min_length=20,
            min_support=2,
            min_mapq=30,
            merge_distance=10,
            threads=10)

        # Run and check result.
        insertions = aligner.extract_insertions(
            read_paths=(Path('/path/to/test.R1.fastq.gz'),
                        Path('/path/to/test.R2.fastq.gz')),
            output_prefix=prefix,
            verbose=False)

        assert len(insertions) == 7

        # Check intermediate calls.
        mock_trim_tr.assert_called_with(
            input_paths=(Path('/path/to/test.R1.fastq.gz'),
                         Path('/path/to/test.R2.fastq.gz')),
            output_paths=(Path(prefix + '.trimmed_tr.R1.fastq.gz'),
                          Path(prefix + '.trimmed_tr.R2.fastq.gz')),
            sequences={'-G': 'TGTATGTAAACTTCCGACTTCAACTG'},
            threads=10,
            verbose=False,
            logger=ANY)

        mock_trim_nt.assert_called_with(
            input_paths=(Path(prefix + '.trimmed_tr.R1.fastq.gz'),
                         Path(prefix + '.trimmed_tr.R2.fastq.gz')),
            output_paths=(Path(prefix + '.trimmed_nt.R1.fastq.gz'),
                          Path(prefix + '.trimmed_nt.R2.fastq.gz')),
            sequences={
                '-a': 'CTGTCTCTTATA',
                '-A': 'CTGTCTCTTATA',
            },
            min_length=20,
            threads=10,
            verbose=False,
            logger=ANY)

        mock_align.assert_called_with(
            read_paths=(Path(prefix + '.trimmed_nt.R1.fastq.gz'),
                        Path(prefix + '.trimmed_nt.R2.fastq.gz')),
            index_path=Path('/path/to/index'),
            output_path=Path(str(subset_bam)),
            threads=10,
            verbose=False,
            logger=ANY)

        mock_extract.assert_called_with(
            alignment_path=Path(str(subset_bam)),
            position_func=ANY,
            sample_func=ANY,
            min_mapq=30,
            paired=True,
            merge_dist=10,
            min_support=2)
