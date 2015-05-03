from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
from builtins import (ascii, bytes, chr, dict, filter, hex, input,
                      int, map, next, oct, open, pow, range, round,
                      str, super, zip)
from future.utils import native_str

import re
from argparse import ArgumentParser
from pathlib import Path

import pandas as pd

from tkgeno.io import GtfFile
from tkgeno.viz.tracks import FeatureTrack, GeneRegionTrack, plot_tracks

try:
    import seaborn as sns
    sns.set_style('whitegrid')
except ImportError:
    sns = None


def setup_parser():
    parser = ArgumentParser(prog='pyim-merge')

    parser.add_argument('input', type=Path)
    parser.add_argument('output', type=Path)
    parser.add_argument('reference_gtf', type=Path)

    parser.add_argument('--region', required=True)
    parser.add_argument('--figure_width', default=None)

    parser.add_argument('--line_height_gene', type=float, default=0.5)
    parser.add_argument('--line_height_insertion', type=float, default=0.5)

    return parser


def main():
    parser = setup_parser()
    args = parser.parse_args()

    # Setup insertion track.
    if sns is not None:
        palette = sns.color_palette()
        color_map = {1: palette[0], -1: palette[2]}
    else:
        color_map = {1: 'red', -1: 'blue'}

    ins_frame = pd.read_csv(str(args.input), sep=native_str('\t'),
                            dtype={'seqname': str})
    ins_track = FeatureTrack.from_location(
        ins_frame, width=2000, line_height=args.line_height_insertion,
        color='strand', color_map=color_map)

    # Setup gene region track.
    gtf_file = GtfFile(args.reference_gtf)
    gene_track = GeneRegionTrack(gtf_file, line_height=args.line_height_gene)

    seqname, start, end = re.split('[:-]', args.region)
    start, end = int(start), int(end)

    fig, axes = plot_tracks([ins_track, gene_track],
                            seqname=seqname, start=start, end=end,
                            figsize=(int(args.figure_width), None))
    fig.savefig(str(args.output), bbox_inches='tight')

if __name__ == '__main__':
    main()
