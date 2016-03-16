from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
from builtins import (ascii, bytes, chr, dict, filter, hex, input,
                      int, map, next, oct, open, pow, range, round,
                      str, super, zip)

import logging
from itertools import chain, repeat

import pandas as pd
from rpy2 import robjects
from rpy2.robjects.packages import importr

from pyim.util.rpy2 import dataframe_to_pandas

from ._util import select_closest


CHROM_MAP = dict(zip(
    list(map(str, range(1, 19+1))) + ['X', 'Y'],
    range(1, 21+1)
))


def register(subparsers, name='kcrbm'):
    parser = subparsers.add_parser(name, help=name + ' help')

    # Required arguments.
    parser.add_argument('input')
    parser.add_argument('output')

    # Optional arguments.
    parser.add_argument('--reference', default='mm10', choices={'mm10'})
    parser.add_argument('--method', default='genes',
                        choices={'genes', 'transcripts'})
    parser.add_argument('--system', default='SB',
                        choices={'MMTV', 'MuLV', 'SB'})
    parser.add_argument('--closest', default=False, action='store_true')

    # Set main for dispatch.
    parser.set_defaults(main=main)

    return parser


def main(args):
    logger = logging.getLogger()

    # Read insertions.
    logger.info('Annotation insertions')
    insertions = pd.read_csv(args.input, sep='\t', dtype={'chrom': str})
    logger.info('Read {} insertions'.format(len(insertions)))

    # Annotate with kcrbm.
    annotation = annotate(insertions, args.reference,
                          args.system, args.method)

    if args.closest:
        # Sub-select for closest features.
        logger.info('Reducing to closest features')
        annotation = select_closest(annotation, col='gene_distance')

    # Merge annotation.
    logger.info('Merging annotation')
    merged = pd.merge(insertions, annotation, on='id', how='left')
    merged.to_csv(args.output, sep='\t', index=False)


def annotate(insertions, reference, system, method):
    # Convert to kcrbm format.
    ins_kcrbm = _convert_to_kcrbm(insertions)

    # Run KCRBM.
    kcrbm = importr('kcrbm')
    genome = _load_genome(reference)

    result = kcrbm.kcrbm(edata=genome, idata=ins_kcrbm, rules=system,
                         reference=reference, map_to=method)
    result = dataframe_to_pandas(result)

    # Convert to gene/transcript frame.
    if method == 'gene':
        result = _convert_gene_result(result)
    elif method == 'transcript':
        result = _convert_transcript_result(result)
    else:
        raise ValueError('Unknown method {}'.format(method))

    return result


def _convert_to_kcrbm(insertion):
    # Extract and rename required columns.
    kcrbm_frame = insertion.ix[:, ['id', 'seqname', 'location', 'strand']]
    kcrbm_frame.columns = ['id', 'chr', 'base', 'ori']

    # Remove any eccentric chromosomes from frame.
    seq_mask = kcrbm_frame.chr.isin(CHROM_MAP.keys())
    if any(~seq_mask):
        dropped_chr = set(kcrbm_frame.ix[~seq_mask].chr)
        print('Warning: dropped insertions not in regular '
              'chromosomes ({})'.format(', '.join(dropped_chr)))

        kcrbm_frame = kcrbm_frame.ix[seq_mask]

    # Convert chr to numeric representation.
    kcrbm_frame['chr'] = kcrbm_frame['chr'].map(CHROM_MAP).astype(int)

    # Copy insertion id to extra column.
    kcrbm_frame['ins_id'] = kcrbm_frame['id']

    return kcrbm_frame


def _load_genome(genome):
    utils = importr("utils")

    if genome == 'mm10':
        utils.data('edata.mm10', package='kcrbm')
        genome_obj = robjects.r['edata.mm10']
    else:
        raise ValueError('Unknown genome version {}'.format(genome))

    return genome_obj


def _convert_gene_result(result):
    result = result.ix[result['ensid'].astype(str) != 'NA']

    gene_distance = result[['d2gss', 'd2gts']]\
        .abs().min(axis=1).astype(int)
    gene_distance.ix[result.mechanism.str.startswith('u')] *= -1

    return pd.DataFrame({
        'insertion_id': result['ins_id'],
        'gene_id': result['ensid'],
        'distance': gene_distance,
        'mechanism': result['mechanism']},
        columns=['insertion_id', 'gene_id', 'distance', 'mechanism'])


def _convert_transcript_result(result):
    result = result.ix[result['ensid'].astype(str) != 'NA']

    transcripts = result['transid'].str.split('|')
    mechanisms = result['mechanism'].str.split('|')

    counts = list(map(len, transcripts))

    ins_id = list(_repeat_list(result['ins_id'], counts))
    ens_id = list(_repeat_list(result['ensid'], counts))

    return pd.DataFrame({'id': ins_id, 'gene': ens_id,
                         'transcript': _flatten_list(transcripts),
                         'mechanism': _flatten_list(mechanisms)},
                        columns=['id', 'gene', 'transcript', 'mechanism'])


def _repeat_list(l, n):
    return chain(*[repeat(el, num) for el, num in zip(l, n)])


def _flatten_list(l):
    return [item for sub_list in l for item in sub_list]
