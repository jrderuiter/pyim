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
    parser.add_argument('--transcript_ids', nargs='+', default=None)

    parser.add_argument('--insertion_width', type=int, default=2000)
    parser.add_argument('--gene_line_height', type=float, default=0.5)
    parser.add_argument('--insertion_line_height', type=float, default=0.5)
    parser.add_argument('--flip_y', default=False, action='store_true')

    parser.add_argument('--dpi', default=72, type=int)

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

    # Parse region.
    seqname, start, end = re.split('[:-]', args.region)
    start, end = int(start), int(end)

    # Setup insertion track.
    ins_frame = pd.read_csv(str(args.input), sep=native_str('\t'),
                            dtype={'seqname': str})
    ins_track = FeatureTrack.from_location(
        ins_frame, width=args.insertion_width,
        line_height=args.insertion_line_height,
        color='strand', color_map=color_map, flip_y=args.flip_y)

    # Setup gene region track.
    gtf = GtfFile(args.reference_gtf)

    if args.transcript_ids is not None:
        gtf = gtf.get_region(seqname, start, end, expand=True)
        gtf = gtf.ix[gtf.transcript_id.isin(args.transcript_ids)]

    gene_track = GeneRegionTrack(gtf, line_height=args.gene_line_height)

    # Setup plot.
    tracks = [gene_track, ins_track] \
        if args.flip_y else [ins_track, gene_track]

    fig, axes = plot_tracks(tracks,
                            seqname=seqname, start=start, end=end,
                            figsize=(int(args.figure_width), None),
                            tick_top=args.flip_y)

    if args.output.suffix == '.png':
        save_kwargs = {'dpi': args.dpi}
    else:
        save_kwargs = {}

    fig.savefig(str(args.output), bbox_inches='tight', **save_kwargs)

if __name__ == '__main__':
    main()
