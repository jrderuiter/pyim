
from __future__ import print_function

import argparse
import subprocess
import pandas as pd
import os, shutil
from os import path
from pyim.insertions import write_to_gff

CHROMOSOMES = ['chr%d' % d for d in range(1,19)] + ['chrX', 'chrY']
CIMPL_SCRIPT = path.join(path.dirname(__file__), '../../r/cimpl.R')


def split_sample_names(sample_names):
    s_list = sample_names.str.split('-').tolist()
    for i in range(len(s_list)):
        s = s_list[i]
        if len(s) == 1:
            s_list[i] = [s[0], 'UN']
        elif len(s) == 2:
            continue
        else: # len > 2
            s_list[i] = [s[0], '-'.join(s[1:])]
    return pd.DataFrame(s_list, columns=['sampleID', 'tissueID'])


def to_cimpl_frame(ins):
    df = split_sample_names(ins['sample'])
    df['tissueType'] = 'Breast'
    df['chr'] = 'chr' + ins['chromosome']
    df['location'] = ins['location']
    df['associatedGene'] = ''
    df['cohort'] = 'sbecad'
    df['contig_depth'] = 1
    df['ins_id'] = ins['ins_id']
    return df


def cimpl(cimpl_frame, workdir, basename):
    file_path = path.abspath(path.join(workdir, basename + '.cimpl-in.tsv'))
    cimpl_frame.to_csv(file_path, sep='\t', index=False)

    output_prefix = path.join(path.abspath(workdir), basename)

    with open(output_prefix + '.cimpl.log', 'w') as log_file:
        subprocess.check_call(['Rscript', CIMPL_SCRIPT, file_path, output_prefix], stderr=log_file)

    cis_list = pd.read_csv(output_prefix + '.cis-list.tsv', sep='\t')
    cis_ins_table = pd.read_csv(output_prefix + '.cis-ins-mapping.tsv', sep='\t')

    return cis_list, cis_ins_table


def merge_cis_list(ins, cis_list, ins_cis_map, prefix='cis_', extra_cols=None):
    if extra_cols is None:
        extra_cols = []

    ins_cis = pd.merge(ins, ins_cis_map, how='right', on='ins_id')

    cis_sel = cis_list[['id', 'pValue', 'scale'] + extra_cols]
    cis_sel.columns = [prefix + c for c in cis_sel.columns]

    return pd.merge(ins_cis, cis_sel, how='left', on='cis_id')


def select_best_cis(ins_cis, cis_prefix='cis_'):
    groupby = ins_cis.groupby('ins_id')

    rows = []
    for key, grp in groupby:
        grp_min_p = grp[grp[cis_prefix + 'pValue'].min() == grp[cis_prefix + 'pValue']]
        grp_min_scale = grp_min_p[grp_min_p[cis_prefix + 'scale'].min() == grp_min_p[cis_prefix + 'scale']]

        if len(grp_min_scale) == 1:
            rows.append(grp_min_scale)
        else:
            print("Warning: skipping insertion with multiple best hits!")

    return pd.concat(rows, ignore_index=True)

def write_cis_bed(path, cis_list):
    with open(path, 'w') as file_:
        print('track name=cis', file=file_)
        cis_list[['chromosome', 'start', 'end', 'id']].to_csv(file_, sep='\t', index=False, header=False)


def main():
    args = setup_parser().parse_args()

    # Setup default arguments
    if args.workdir is None:
        args.workdir = path.join(args.outputdir, '_pyim')
    if not path.exists(args.workdir): os.makedirs(args.workdir)

    if args.basename is None:
        args.basename = path.basename(args.input).split('.')[0]

    # Read and convert frame to CIMPLs format
    ins_frame = pd.read_csv(args.input, sep='\t')
    ins_frame = ins_frame[['ins_id', 'chromosome', 'location', 'lp', 'unique_lp', 'strand', 'sample']]

    cimpl_frame = to_cimpl_frame(ins_frame)
    cimpl_frame = cimpl_frame[cimpl_frame['chr'].isin(CHROMOSOMES)]

    # Run cimpl
    cis_list, cis_ins_table = cimpl(cimpl_frame, args.workdir, args.basename)

    # Merge into insertion table
    extra_cis_cols = ['associatedEnsemblGeneId', 'otherEnsemblGeneId',
                      'associatedExternalGeneId', 'otherExternalGeneId']
    merged_ins = merge_cis_list(ins_frame, cis_list, cis_ins_table, prefix='cis_', extra_cols=extra_cis_cols)

    # Select 'best' CISs
    merged_ins_best = select_best_cis(merged_ins)

    # Write tables
    output_base = path.join(args.outputdir, args.basename)

    cis_list.to_csv(output_base + '.cis.tsv', sep='\t', index=False, float_format='%3.20f')
    merged_ins.to_csv(output_base + '.cis-ins.tsv', sep='\t', index=False, float_format='%3.20f')
    merged_ins_best.to_csv(output_base + '.cis-ins.best.tsv', sep='\t', index=False, float_format='%3.20f')
    write_to_gff(merged_ins_best, output_base + '.cis-ins.best.gff')

    # Write bed files
    cis_best = cis_list[cis_list['id'].isin(merged_ins_best['cis_id'])]
    write_cis_bed(output_base + '.cis.bed', cis_list)
    write_cis_bed(output_base + '.cis-best.bed', cis_best)

    # Clean up our work directory
    #shutil.rmtree(args.workdir, ignore_errors=True)


def setup_parser(root=None):
    if root is not None:
        parser = root.add_parser('cimpl')
    else:
        parser = argparse.ArgumentParser()

    parser.add_argument('-i', '--input', required=True)
    parser.add_argument('-o', '--outputdir', default='.')
    parser.add_argument('-d', '--workdir', default=None)
    parser.add_argument('-b', '--basename', default=None)

    return parser


if __name__ == '__main__':
    main()