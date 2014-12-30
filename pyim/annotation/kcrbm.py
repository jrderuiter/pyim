
import readline  # Work-around for: ../lib/libreadline.so.6: undefined symbol: PC

import pandas
import pandas.rpy.common as com
from numpy import issubdtype
from rpy2 import robjects
from rpy2.robjects.packages import importr

from pyim_common.util.list import repeat_list, flatten_list

from pyim.annotation.annotator import Annotator


class KcrbmAnnotator(Annotator):

    def __init__(self, genome, system):
        super(KcrbmAnnotator, self).__init__()
        if system not in {'SB'}:
            raise ValueError('Unknown system {}'.format(system))

        if genome not in {'mm10'}:
            raise ValueError('Unsupported genome {}'.format(genome))

        self.genome = genome
        self.system = system

    def annotate_by_gene(self, insertions):
        kcrbm_ins = self._convert_to_kcrbm_frame(insertions)
        kcrbm_result = self._run_kcrbm(kcrbm_ins, method='genes')
        return self._parse_gene_result(kcrbm_result)

    def annotate_by_transcript(self, insertions):
        kcrbm_ins = self._convert_to_kcrbm_frame(insertions)
        kcrbm_result = self._run_kcrbm(kcrbm_ins, method='transcripts')
        return self._parse_transcript_result(kcrbm_result)

    def _run_kcrbm(self, kcrbm_frame, method):
        kcrbm_frame['ins_id'] = kcrbm_frame['id']
        ins_kcrbm_r = com.convert_to_r_dataframe(kcrbm_frame)

        kcrbm = importr('kcrbm')
        genome = self._load_genome(self.genome)
        res = kcrbm.kcrbm(edata=genome, idata=ins_kcrbm_r,
                          rules=self.system, reference=self.genome, map_to=method)

        return com.convert_robj(res)

    @staticmethod
    def _convert_to_kcrbm_frame(insertions):
        # Extract and rename required columns.
        kcrbm_frame = insertions.ix[:, ['id', 'seqname', 'location', 'strand']]
        kcrbm_frame.columns = ['id', 'chr', 'base', 'ori']

        # Map chromosomes to numeric values.
        chr_map = dict(zip(
            list(map(str, range(1, 19 + 1))) + ['X', 'Y'],
            range(1, 21 + 1)
        ))
        kcrbm_frame['chr'] = kcrbm_frame['chr'].map(chr_map).astype(int)

        # Convert orientation if required.
        if not issubdtype(kcrbm_frame['ori'].dtype, int):
            kcrbm_frame['ori'] = kcrbm_frame['ori'].map({'+': 1, '-': -1})

        return kcrbm_frame

    @staticmethod
    def _load_genome(genome):
        utils = importr("utils")

        if genome == 'mm10':
            utils.data('edata.mm10', package='kcrbm')
            genome_obj = robjects.r['edata.mm10']
        else:
            raise ValueError('Unknown genome version {}'.format(genome))

        return genome_obj

    @staticmethod
    def _parse_gene_result(result):
        result = result.ix[result['ensid'].astype(str) != 'NA']
        return pandas.DataFrame({'id': result['ins_id'],
                                 'gene': result['ensid'],
                                 'mechanism': result['mechanism']},
                                columns=['id', 'gene', 'mechanism'])

    @staticmethod
    def _parse_transcript_result(result):
        result = result.ix[result['ensid'].astype(str) != 'NA']

        tr_list = result['transid'].str.split('|')
        mech_list = result['mechanism'].str.split('|')

        counts = list(map(len, tr_list))

        ins_id = list(repeat_list(result['ins_id'], counts))
        ens_id = list(repeat_list(result['ensid'], counts))

        return pandas.DataFrame({'id': ins_id, 'gene': ens_id,
                                 'transcript': flatten_list(tr_list),
                                 'mechanism': flatten_list(mech_list)},
                                columns=['id', 'gene', 'transcript', 'mechanism'])
