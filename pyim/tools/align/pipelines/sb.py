import pandas
from skbio import DNASequence
from skbio.io import read as skbio_read
from pyim_common.model.insertion import Insertion

from pyim.common.alignment.vector.aligners import VectorAligner
from .base import Pipeline
from ..pyim.tools.alignment.genome.aligners import Bowtie2Aligner


class SbPipeline(Pipeline):

    def __init__(self, reference, barcode_file, sb_sequence,
                 t7_sequence, min_genomic_length=15):
        super(SbPipeline, self).__init__()

        self.reference = reference
        self.barcodes = list(skbio_read(barcode_file, format='fasta',
                                        constructor=DNASequence))[:10]

        self.sb_sequence = DNASequence(sb_sequence, id='SB')
        self.t7_sequence = DNASequence(t7_sequence, id='T7')

        self.min_genomic_length = min_genomic_length

    @classmethod
    def configure_argparser(cls, parser):
        parser = super(SbPipeline, cls).configure_argparser(parser)

        parser.add_argument('reference')
        parser.add_argument('barcode_file')

        parser.add_argument('--sb_sequence', default='GTGTATGTAAACTTCCGACTTCAAC')
        parser.add_argument('--t7_sequence', default='CCTATAGTGAGTCGTATTA')
        parser.add_argument('--min_genomic_length', default=15, type=int)

        return parser

    def run(self, input_path, output_path):
        sequences = skbio_read(input_path, format='fasta', constructor=DNASequence)

        # Match barcodes
        barcodes = match_barcodes(sequences, self.barcodes)

        # Align to vectors
        # sb_aligner = SequenceAlignerFactory.make_from_options(config['sb_aligner'])

        # t7 = DNASequence(config['t7_sequence'], id='T7')
        # t7_aligner = SequenceAlignerFactory.make_from_options(config['t7_aligner'])

        # vec_alignments = align_to_vectors(sequences, [sb, t7], [sb_aligner, t7_aligner])

        # Extract genomic sequences using vector alignments.
        genomic_seqs = extract_genomic(sequences, vec_alignments)
        genomic_seqs = {k: v for k, v in genomic_seqs.items()
                        if len(v.sequence) >= config['min_genomic_length']}

        genomic_alns = align_to_reference(genomic_seqs.values(), config['reference'])

        # Map alignments to unique insertions.
        insertions = map_insertions(genomic_alns, barcodes)
        Insertion.to_file(insertions, output_path)


def match_barcodes(sequences, barcodes):
    bc_aligner = VectorAligner()
    alignments = bc_aligner.align_multiple(sequences, barcodes)
    return {aln.q_name: aln.v_name for aln in alignments.values()}


def align_to_vectors(sequences, vectors, aligners=None):
    if aligners is None:
        aligners = [VectorAligner() for _ in vectors]

    alignments = {}
    for vector, aligner in zip(vectors, aligners):
        alignments[vector.name] = aligner.align(sequences, vector)
    
    return alignments


def extract_genomic(sequences, alignments):
    sb_alns, t7_alns = alignments['SB'], alignments['T7']
    
    genomic_seqs = {}
    for seq in sequences:
        if seq.name in sb_alns and seq.name in t7_alns:
            start = sb_alns[seq.name].q_end
            end = t7_alns[seq.name].q_start

            new_seq = DNASequence(seq.sequence[start:end], id=seq.name)

            genomic_seqs[seq.name] = new_seq

    return genomic_seqs


def align_to_reference(sequences, reference):
    aligner = Bowtie2Aligner(reference)
    alignments = aligner.align(sequences, work_dir='/Users/Julian/Desktop/Bowtie')
    return alignments


def map_insertions(genomic_alignments, barcodes):
    # Convert alignments to a data frame for convenient aggregation.
    aln_frame = pandas.DataFrame((ga._asdict() for ga in genomic_alignments.values()))
    aln_frame['barcode'] = aln_frame['q_name'].map(barcodes)
    aln_frame.dropna(subset=['barcode'], inplace=True)

    # Aggregate insertions per strand.
    fwd_ins = _group_per_strand(aln_frame, 1, 'r_start', 'r_end')
    rev_ins = _group_per_strand(aln_frame, -1, 'r_end', 'r_start', id_offset=len(fwd_ins))

    return fwd_ins + rev_ins


def _group_per_strand(aln_frame, strand, anchor_field, var_field, id_offset=0):
    alns = aln_frame.ix[aln_frame['r_strand'] == strand]
    
    # Count the number of (unique) T7 ligation points
    agg_funcs = {'lp': len, 'unique_lp': lambda x: len(x.unique())}

    grouped = alns.groupby(['barcode', 'r_name', anchor_field])
    group_stats = grouped[var_field].agg(agg_funcs)

    # Create insertion frame from statistics
    group_stats = group_stats.reset_index()

    insertions = []
    for i, (_, row) in enumerate(group_stats.iterrows()):
        name = 'INS.{}'.format((i + id_offset) + 1)
        metadata = {'lp': row.lp, 'uniqueLp': row.unique_lp}
        ins = Insertion(name=name, seqname=row.r_name, location=row[anchor_field],
                        strand=strand, barcode=row.barcode, metadata=metadata)
        insertions.append(ins)
    
    return insertions
