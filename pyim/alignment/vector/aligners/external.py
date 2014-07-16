
import tempfile
import uuid
import shutil
import subprocess
from os import path

from pyim.io import write_fasta
from pyim.util import chunks
from pyim.alignment.vector.aligners.base import InexactReadAligner
from pyim.alignment.vector.model import Alignment


class ExonerateReadAligner(InexactReadAligner):

    def __init__(self, min_score=0, model='affine:local', exhaustive=False, bestn=1, filters=None):
        super(ExonerateReadAligner, self).__init__(filters=filters, min_score=min_score)
        self.model = model
        self.exhaustive = exhaustive
        self.bestn = bestn

    def _setup(self):
        tmp_dir = tempfile.mkdtemp()

        path_template = path.join(tmp_dir, '{uuid}.{ext}')
        read_path = path_template.format(uuid=uuid.uuid4(), ext='fa')
        target_path = path_template.format(uuid=uuid.uuid4(), ext='fa')
        output_path = path_template.format(uuid=uuid.uuid4(), ext='exo')

        return read_path, target_path, output_path, tmp_dir

    def _cleanup(self, tmp_dir):
        try: shutil.rmtree(tmp_dir)
        except OSError as exc:
            if exc.errno != 2: raise

    def _align_target(self, reads, target):
        # Setup and write sequences
        read_path, target_path, output_path, tmp_dir = self._setup()
        write_fasta(read_path, reads)
        write_fasta(target_path, [target])

        # Do your thing, parse result and clean up
        self._run(read_path, target_path, output_path)
        alignments, unmapped = self._parse_result(output_path, reads, target)
        self._cleanup(tmp_dir)

        return alignments, unmapped

    def _run(self, read_path, target_path, output_path):
        exo_cmd = "exonerate -s {minscore} -m {model} --showalignment 0 --showvulgar 0 --bestn {bestn} " + \
                            "--ryo '%qi %qab %qae %ti %tab %tae %s %pi %qS %V\n' {reads} {target} > {output}"
        if self.exhaustive: exo_cmd += ' -E'

        exo_cmd_frmt = exo_cmd.format(minscore=self.min_score, model=self.model, bestn=self.bestn,
                                      reads=read_path, target=target_path, output=output_path)

        with open('/dev/null', 'w') as dev_null:
            subprocess.check_call(exo_cmd_frmt, shell=True, stderr=dev_null)

    def _parse_result(self, output_path, reads, target):
        with open(output_path, 'r') as file_:
            lines = file_.readlines()[3:-1]
            split_lines = [l.strip().split(' ') for l in lines]

        read_lookup = { r.name: r for r in reads }
        alignment_type = 'exonerate[%d' % self.min_score
        alignment_type += ']' if not self.exhaustive else ',exhaustive]'

        alignments = []
        for sl in split_lines:
            if sl[8] == '+':
                q_start, q_seq = int(sl[1]), read_lookup[sl[0]].seq
                t_start, t_seq = int(sl[4]), target.seq
                identity = float(sl[7])

                aln_str = self._parse_vulgar(sl[9:], identity, ((q_start, q_seq), (t_start, t_seq)))

                alignments.append(Alignment(query_name=sl[0],  query_start=q_start,  query_end=int(sl[2]),
                                            target_name=sl[3], target_start=t_start, target_end=int(sl[5]),
                                            score=int(sl[6]), identity=identity, type=alignment_type,
                                            alignment=aln_str, query_seq=q_seq, target_seq=t_seq))

        alignments = self._alignments_to_frame(alignments)
        unmapped_ids = set(read_lookup.keys()).difference(set(alignments['query_name']))
        unmapped = [read_lookup[i] for i in unmapped_ids]

        return alignments, unmapped

    def _parse_vulgar(self, vulgar, identity, alignment_info):
        vulgar_blocks = [tuple([c[0], int(c[1]), int(c[2])]) for c in chunks(vulgar, 3)]

        if identity < 100.0:
            aln_str = self._parse_vulgar_mismatched(vulgar_blocks, alignment_info)
        else:
            aln_str = self._parse_vulgar_matched(vulgar_blocks)

        return aln_str

    def _parse_vulgar_matched(self, vulgar_blocks):
        alignment_str = ''
        for aln_type, query_len, target_len in vulgar_blocks:
            if aln_type == 'M':
                    alignment_str += query_len * 'M'
            elif aln_type == 'G':
                if query_len > 0:
                    alignment_str += query_len * 'D'
                elif target_len > 0:
                    alignment_str += target_len * 'I'
            else: raise ValueError('Unknown type in vulgar block')
        return alignment_str

    def _parse_vulgar_mismatched(self, vulgar_blocks, alignment_info):
        (q_start, q_seq), (t_start, t_seq) = alignment_info
        qi, ti = q_start, t_start

        alignment_str = ''
        for aln_type, query_len, target_len in vulgar_blocks:
            if aln_type == 'M':
                for i in range(query_len):
                    alignment_str += 'M' if q_seq[qi+i] == t_seq[ti+i] else 'S'
                qi, ti = qi + query_len, ti + query_len
            elif aln_type == 'G':
                if query_len > 0:
                    alignment_str += query_len * 'D'
                    qi += query_len
                elif target_len > 0:
                    alignment_str += query_len * 'I'
                    ti += target_len
            else: raise ValueError('Unknown type in vulgar block')

        return alignment_str