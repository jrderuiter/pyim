from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
from builtins import (ascii, bytes, chr, dict, filter, hex, input,
                      int, map, next, oct, open, pow, range, round,
                      str, super, zip)

from itertools import chain, repeat
from pathlib import Path

import pandas as pd
from numpy import issubdtype
from rpy2 import robjects

from tkgeno.util.rpy2 import importr, pandas_to_dataframe, dataframe_to_pandas

from .base import Annotator


class KcRbmAnnotator(Annotator):

    def __init__(self, reference, system):
        super().__init__()

        if system not in {'SB'}:
            raise ValueError('Unknown system {}'.format(system))

        if reference not in {'mm10'}:
            raise ValueError('Unsupported genome {}'.format(reference))

        self._reference = reference
        self._system = system

    @classmethod
    def configure_argparser(cls, subparsers, name='kcrbm'):
        parser = subparsers.add_parser(name, help=name + ' help')

        parser.add_argument('input', type=Path)
        parser.add_argument('output', type=Path)

        parser.add_argument('--reference', default='mm10')
        parser.add_argument('--system', default='SB')

        return parser

    @classmethod
    def from_args(cls, args):
        return cls(reference=args['reference'], system=args['system'])

    def annotate(self, frame, type_='gene'):
        kcrbm_ins = self._convert_to_kcrbm_frame(frame)
        kcrbm_result = self._run_kcrbm(kcrbm_ins, method='genes')

        gene_mapping = self._parse_gene_result(kcrbm_result)
        return pd.merge(frame, gene_mapping, on='insertion_id')

    @staticmethod
    def _convert_to_kcrbm_frame(frame):
        # Extract and rename required columns.
        kcrbm_frame = frame.ix[:, ['insertion_id', 'seqname',
                                   'location', 'strand']]
        kcrbm_frame.columns = ['id', 'chr', 'base', 'ori']

        # Map chromosomes to numeric values.
        chr_map = dict(zip(
            list(map(str, range(1, 19+1))) + ['X', 'Y'],
            range(1, 21+1)
        ))
        kcrbm_frame['chr'] = kcrbm_frame['chr'].map(chr_map).astype(int)

        # Convert orientation if required.
        if not issubdtype(kcrbm_frame['ori'].dtype, int):
            kcrbm_frame['ori'] = kcrbm_frame['ori'].map({'+': 1, '-': -1})

        return kcrbm_frame

    def _run_kcrbm(self, kcrbm_frame, method):
        kcrbm_frame['ins_id'] = kcrbm_frame['id']
        kcrbm_df = pandas_to_dataframe(kcrbm_frame)

        kcrbm = importr('kcrbm')
        genome = self._load_genome(self._reference)
        res = kcrbm.kcrbm(edata=genome, idata=kcrbm_df, rules=self._system,
                          reference=self._reference, map_to=method)

        return dataframe_to_pandas(res)

    @staticmethod
    def _parse_gene_result(result):
        result = result.ix[result['ensid'].astype(str) != 'NA']
        return pd.DataFrame({'insertion_id': result['ins_id'],
                             'gene_id': result['ensid'],
                             'mechanism': result['mechanism']},
                            columns=['insertion_id', 'gene_id', 'mechanism'])

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
    def _parse_transcript_result(result):
        result = result.ix[result['ensid'].astype(str) != 'NA']

        tr_list = result['transid'].str.split('|')
        mech_list = result['mechanism'].str.split('|')

        counts = list(map(len, tr_list))

        ins_id = list(_repeat_list(result['ins_id'], counts))
        ens_id = list(_repeat_list(result['ensid'], counts))

        return pd.DataFrame({'id': ins_id, 'gene': ens_id,
                             'transcript': _flatten_list(tr_list),
                             'mechanism': _flatten_list(mech_list)},
                            columns=['id', 'gene', 'transcript', 'mechanism'])


def _repeat_list(l, n):
    return chain(*[repeat(el, num) for el, num in zip(l, n)])


def _flatten_list(l):
    return [item for sub_list in l for item in sub_list]