
import subprocess, pandas

from os import path
from pyim.genomic.base import GenomicAligner, GenomicAlignment

PSL_COLUMN_NAMES = ['match', 'mismatch', 'repmatch', 'num_n',
                    'q_gap_count', 'q_gap_bases',
                    't_gap_count', 't_gap_bases', 'strand',
                    'q_name', 'q_size', 'q_start', 'q_end',
                    't_name', 't_size', 't_start', 't_end',
                    'block_count', 'block_sizes',
                    'q_starts', 't_starts']


def read_psl(psl_file):
    return pandas.read_table(psl_file, header=None, skiprows=5, names=PSL_COLUMN_NAMES)


class BlatAligner(GenomicAligner):

    def __init__(self, filters=None, minscore=10):
        super(BlatAligner, self).__init__(filters)
        self.minscore = minscore
        
    def _run(self, reads, reference, tmp_dir):
        read_file = self._write_fasta(tmp_dir, reads)

        output_path = path.join(tmp_dir, 'mapped_reads.psl')
        log_path = path.join(tmp_dir, 'blat.log')

        cmd = "blat {reference} {fasta} -out=psl -minScore={minscore} {psl} 2> {log}"
        cmd_frmt = cmd.format(reference=reference, fasta=read_file,
                              minscore=self.minscore, psl=output_path, log=log_path)
        subprocess.check_call(cmd_frmt, shell=True)

        return output_path

    def _read_output(self, psl_file):
        psl_frame = read_psl(psl_file)
        alignments = [GenomicAlignment(query_name=row['q_name'],
                                       query_start=row['q_start'],
                                       query_end=row['q_end'],
                                       target_name=row['t_name'],
                                       target_start=row['t_start'],
                                       target_end=row['t_end'],
                                       strand=row['strand'],
                                       alignment='') for row in psl_frame.iterrows()]

        return alignments
