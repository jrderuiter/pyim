import pandas

from pyim.common.model import Sequence
from pyim.alignment.general.vector.aligners import VectorAligner
from pyim.alignment.genome.aligners import Bowtie2Aligner
from pyim.common.io import read_fasta
from pyim.common.model import Insertion

SB_SEQUENCE = 'GTGTATGTAAACTTCCGACTTCAAC'
T7_SEQUENCE = 'CCTATAGTGAGTCGTATTA'

REFERENCE = '/Volumes/Datastore/Julian/References/mus_musculus/' + \
            'mm10/dna/Bowtie2/Mus_musculus.GRCm38.dna.primary_assembly'

BARCODES = list(read_fasta('/Users/Julian/Software/python/im/pyim_new/data/SB_barcodes.fa'))[0:10]

class SbPipeline(object):

    @classmethod
    def run(Class, options):
        # TODO: supply options externally (and use them).
        # TODO: somehow specify which aligners are used (and their options).

        sequences = list(read_fasta(options.input))

        ## Match barcodes
        barcodes = match_barcodes(sequences, BARCODES)

        ## Align to vectors
        sb = Sequence(name='SB', sequence=SB_SEQUENCE)
        t7 = Sequence(name='T7', sequence=T7_SEQUENCE)

        vec_alignments = align_to_vectors(sequences, [sb, t7])

        ## Extract genomic sequences using vector alignments.
        genomic_seqs = extract_genomic(sequences, vec_alignments)
        genomic_seqs = {k: v for k, v in genomic_seqs.items() if len(v.sequence) >= 15}

        genomic_alns = align_to_reference(genomic_seqs.values(), REFERENCE)

        ## Map alignments to unique insertions.
        insertions = map_insertions(genomic_alns, barcodes)
        Insertion.to_file(insertions, '/Users/Julian/Desktop/insertions.txt')


def match_barcodes(sequences, barcodes):
    bc_aligner = VectorAligner()
    alignments = bc_aligner.align_multiple(sequences, barcodes)
    return {aln.q_name: aln.v_name for aln in alignments.values()}


def align_to_vectors(sequences, vectors, aligners=None):
    if aligners is None:
        aligners = [VectorAligner() for v in vectors]

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

            new_seq = Sequence(name=seq.name, 
                               sequence=seq.sequence[start:end])

            genomic_seqs[seq.name] = new_seq

    return genomic_seqs


def align_to_reference(sequences, reference):
    aligner = Bowtie2Aligner(reference)
    alignments = aligner.align(sequences, work_dir='/Users/Julian/Desktop/Bowtie')
    return alignments


def map_insertions(genomic_alignments, barcodes):
    ## Convert alignments to a data frame for convenient aggregation.
    aln_frame = pandas.DataFrame((ga._asdict() for ga in genomic_alignments.values()))
    aln_frame['barcode'] = aln_frame['q_name'].map(barcodes)
    aln_frame.dropna(subset=['barcode'], inplace=True)

    ## Aggregate insertions per strand.
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
