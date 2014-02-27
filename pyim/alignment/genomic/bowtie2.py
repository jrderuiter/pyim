import subprocess

import pysam

from pyim.alignment.genomic.base import GenomicAligner, GenomicAlignment


class Bowtie2Aligner(GenomicAligner):

    def __init__(self, filters=None, max_hits=2, num_cores=1):
        super(Bowtie2Aligner, self).__init__(filters)
        self.num_cores = num_cores
        self.max_hits = max_hits


    def _run(self, reads, reference, tmp_dir):
        read_file = self._write_fasta(tmp_dir, reads)
        output_path = self._tmpfile(tmp_dir, ext='sam')

        cmd = "bowtie2 -p {n_cpus} -k {n_hits} -x {reference} -f -U {fasta} -S {output}"
        cmd_frmt = cmd.format(reference=reference, fasta=read_file, n_hits=self.max_hits,
                              output=output_path, n_cpus=self.num_cores)
        subprocess.check_call(cmd_frmt, shell=True)

        return output_path

    def _read_output(self, file_path):
        sam_file = pysam.Samfile(file_path, "r" )

        alignments= []
        for read in sam_file.fetch():
            if read.tid != -1:
                alignment = GenomicAlignment(query_name=read.qname,
                                             query_start=read.qstart,
                                             query_end=read.qend,
                                             query_size=read.rlen,
                                             target_name=sam_file.getrname(read.tid),
                                             target_start=read.positions[0],
                                             target_end=read.positions[-1],
                                             strand='-' if read.is_reverse else '+',
                                             alignment=read.cigarstring)
                alignments.append(alignment)

        return alignments
