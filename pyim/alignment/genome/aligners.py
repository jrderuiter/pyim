import shutil
import subprocess
import tempfile
from os import path

import pysam

from pyim.alignment.genome.model import GenomicAlignment
from pyim.common.io import makedirs_safe, write_fasta


class ReferenceAligner(object):

    def __init__(self, reference, options=None):
        if options is None:
            options = {}
            
        self.reference = reference
        self.options = options

    def align(self, sequences, work_dir=None, keep_work=False):
        work_dir = self._setup_work_dir(work_dir)
        seq_path = self._write_sequences(sequences, work_dir)
        
        aln_path = self._run(seq_path, work_dir)
        alignments = self._parse_alignments(aln_path)

        if not keep_work:
            self._cleanup_work_dir(work_dir)

        return alignments

    def _run(self, seq_path, work_dir):
        raise NotImplementedError

    def _parse_alignments(self, aln_path):
        raise NotImplementedError

    def _write_sequences(self, sequences, work_dir, 
                         file_name='sequences.fna'):
        seq_path = path.join(work_dir, file_name)
        write_fasta(sequences, seq_path)
        return seq_path

    def _setup_work_dir(self, work_dir):
        if work_dir is None:
            work_dir = tempfile.mkdtemp()
        else:
            makedirs_safe(work_dir)
        return work_dir

    def _cleanup_work_dir(self, work_dir):
        shutil.rmtree(work_dir)


class Bowtie2Aligner(ReferenceAligner):

    def _run(self, seq_path, work_dir):
        aln_path = path.join(work_dir, 'alignments.sam')

        ## Setup required arguments.
        args = ['-x', self.reference, 
                '-f', '-U', seq_path, 
                '-S', aln_path,
                '-k', '1']

        ## Add any optional arguments if present.
        if 'threads' in self.options:
            args += ['-p', str(self.options['threads'])]

        ## Actually run bowtie2!
        subprocess.check_call(['bowtie2'] + args)

        return aln_path

    def _parse_alignments(self, aln_path):
        sam_file = pysam.Samfile(aln_path, 'r')

        alignments = {}
        for read in sam_file.fetch():
            if read.tid != -1:
                q_name = read.qname

                if q_name in alignments:
                    raise ValueError('Encountered multiple alignments ' + \
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
