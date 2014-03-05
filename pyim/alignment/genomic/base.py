
import tempfile
import shutil
import uuid
import logging
import pandas

from os import path
from collections import namedtuple

from pyim.io import write_fasta


GenomicAlignment = namedtuple('GenomicAlignment',
                              ['query_name', 'query_start', 'query_end', 'query_size',
                               'target_name', 'target_start', 'target_end',
                               'strand', 'alignment'])


class GenomicAligner(object):

    def __init__(self, filters=None):
        if filters is None: filters = []
        self.filters = filters
        self.logger = logging.getLogger(self.__class__.__name__)

    def align(self, reads, reference):
        tmp_dir = self._setup()

        output_path = self._run(reads, reference, tmp_dir)
        alignments = self._read_output(output_path)
        alignment_df = pandas.DataFrame(alignments, columns=alignments[0]._fields)

        filtered_alignment = self.filter(alignment_df)
        dupl, uniq = self.duplicates(filtered_alignment)
        unaligned = self.unaligned(reads, uniq, dupl)

        self._cleanup(tmp_dir)
        self._stats(reads, uniq, dupl, unaligned)

        hits = self.to_hits(uniq)
        return hits, uniq, dupl, unaligned

    def _run(self, reads, reference, tmp_dir):
        raise NotImplementedError

    def _read_output(self, output_path):
        raise NotImplementedError

    def filter(self, alignment):
        filtered_alignment = alignment
        for filt in self.filters:
            filtered_alignment, _ = filt.apply(filtered_alignment)
        return filtered_alignment

    def duplicates(self, alignments):
        duplicates = alignments['query_name'][alignments.duplicated('query_name')]
        duplicate_mask = ~alignments['query_name'].isin(duplicates)

        dupl_alignments = alignments[~duplicate_mask]
        uniq_alignments = alignments[duplicate_mask]

        return dupl_alignments, uniq_alignments

    def unaligned(self, reads, uniq_alignments, dupl_alignments):
        read_names = pandas.Series([r.name for r in reads])
        read_unaligned = ~read_names.isin(uniq_alignments['query_name']) & \
                         ~read_names.isin(dupl_alignments['query_name'])
        unaligned_reads = [reads[i] for i, unaligned in enumerate(read_unaligned) if unaligned]
        return unaligned_reads

    def _stats(self, reads, uniq, dupl, unaligned):
        n_uniq, n_dupl = len(uniq), dupl['query_name'].nunique()
        n_reads, n_unaligned = float(len(reads)), len(unaligned)

        self.logger.info('%d reads (%02.2f%%) aligned exactly 1 time' % (n_uniq, (n_uniq/n_reads)*100))
        self.logger.info('%d reads (%02.2f%%) aligned >1 times' % (n_dupl, (n_dupl/n_reads)*100))
        self.logger.info('%d reads (%02.2f%%) aligned 0 times' % (n_unaligned, (n_unaligned/n_reads)*100))

    def _setup(self):
        tmp_dir = tempfile.mkdtemp()
        return tmp_dir

    def _tmpfile(self, tmp_dir, ext):
        fn = '{fn}.{ext}'.format(fn=uuid.uuid4(), ext=ext)
        return path.join(tmp_dir, fn)

    def _cleanup(self, tmp_dir):
        try: shutil.rmtree(tmp_dir)
        except OSError as exc:
            if exc.errno != 2: raise

    def _write_fasta(self, dir_path, reads):
        file_path = self._tmpfile(dir_path, ext='fna')
        return write_fasta(file_path, reads)

    def to_hits(self, alignment_df):
        transformed = alignment_df[['query_name', 'target_name', 'target_start', 'target_end', 'strand']]
        transformed.columns = ['query_name', 'contig_name', 'start', 'end', 'strand']
        return transformed

 

