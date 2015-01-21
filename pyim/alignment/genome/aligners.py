import subprocess
from os import path

import pysam

from skbio.io import write as skbio_write

from pyim_common.io import work_directory

from .model import GenomicAlignment


class ReferenceAligner(object):

    def __init__(self, reference, options=None):
        if options is None:
            options = {}
            
        self.reference = reference
        self.options = options

    def align(self, sequences, work_dir=None, keep_work=False):

        with work_directory(work_dir, keep=keep_work) as w_dir:
            seq_path = self._write_sequences(sequences, work_dir)
        
            aln_path = self._run(seq_path, w_dir)
            alignments = self._parse_alignments(aln_path)

        return alignments

    def _run(self, seq_path, work_dir):
        raise NotImplementedError

    def _parse_alignments(self, aln_path):
        raise NotImplementedError

    @classmethod
    def _write_sequences(cls, sequences, work_dir, file_name='sequences.fna'):

        seq_path = path.join(work_dir, file_name)
        skbio_write(sequences, format='fasta', into=seq_path)

        return seq_path


class Bowtie2Aligner(ReferenceAligner):

    def _run(self, seq_path, work_dir):
        aln_path = path.join(work_dir, 'alignments.sam')

        # Setup required arguments.
        args = ['-x', self.reference, 
                '-f', '-U', seq_path, 
                '-S', aln_path,
                '-k', '1']

        # Add any optional arguments if present.
        if 'threads' in self.options:
            args += ['-p', str(self.options['threads'])]

        # Actually run bowtie2!
        subprocess.check_call(['bowtie2'] + args)

        return aln_path

    def _parse_alignments(self, aln_path):
        sam_file = pysam.Samfile(aln_path, 'r')

        alignments = {}
        for read in sam_file.fetch():
            if read.tid != -1:
                q_name = read.qname

                if q_name in alignments:
                    raise ValueError('Encountered multiple alignments '
                                     'for read %s' % q_name)

                aln = GenomicAlignment(
                    q_name=q_name, 
                    q_start=read.qstart, 
                    q_end=read.qend,
                    q_length=read.rlen,
                    r_name=sam_file.getrname(read.tid), 
                    r_start=read.positions[0], 
                    r_end=read.positions[-1], 
                    r_strand=-1 if read.is_reverse else 1)

                alignments[q_name] = aln

        return alignments
