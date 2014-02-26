
import pandas, logging

from pyim.alignment.model import Alignment
from pyim.io import Sequence
from pyim.vector.base import VectorAligner, AlignmentResult
from pyim.ggplot import ggplot_bar_hist, ggplot_pie, ggplot_stacked_bar, ggplot_hist


class SBVectorAligner(VectorAligner):

    def __init__(self):
        super(self.__class__, self).__init__()
        self.logger = logging.getLogger('VectorAligner')

    def align(self, reads, vector_file, barcode_file,
                    sb_aligner=None, t7_aligner=None, bc_aligner=None, missing_T7=False):
        vectors = self._read_seqs(vector_file, name_mask=['T7', 'SB'])
        barcodes = self._read_seqs(barcode_file)

        self.logger.info('Aligning SB and T7 vectors')
        sb_alignments, sb_unmapped = self._align_sb_vector(reads, vectors['SB'], sb_aligner)
        t7_alignments, t7_unmapped = self._align_t7_vector(reads, vectors['T7'], t7_aligner, missing_T7)

        self.logger.info('Aligning barcode vectors')
        bc_alignments, bc_unmapped = self._align_barcodes(reads, barcodes, bc_aligner)

        self.logger.info('Combining alignments')
        combined_alignment = self._combine_alignments(sb_alignments, t7_alignments, bc_alignments)
        combined_alignment.to_csv('vector_alignments.tsv', sep='\t')

        self.logger.info('Extracting genomic sequences')
        gen_sequences = self._genomic_sequences(combined_alignment)

        result = AlignmentResult(combined_alignments=combined_alignment, genomic_sequences=gen_sequences,
                                 vector_alignments={ 'T7': t7_alignments, 'SB': sb_alignments },
                                 vector_unmapped={ 'T7': t7_unmapped, 'SB': sb_unmapped },
                                 barcode_alignments=bc_alignments, barcode_unmapped=bc_unmapped)
        return result

    def _align_sb_vector(self, reads, sb_vector, aligner):
        return self._align_vector(reads, sb_vector, aligner)

    def _align_t7_vector(self, reads, t7_vector, aligner, allow_missing_t7=False):
        t7_alignments, t7_unmapped = self._align_vector(reads, t7_vector, aligner)

        if allow_missing_t7:
            dummy_alns = [Alignment(r.name, len(r.seq), len(r.seq), r.seq,
                                    'T7', 0, 0, t7_vector.seq, 0, 0, '', 'dummy') for r in t7_unmapped]
            dummy_aln_frame = _alignments_to_frame(dummy_alns)
            t7_alignments = pandas.concat([t7_alignments, dummy_aln_frame], ignore_index=True)
            t7_unmapped = []

        return t7_alignments, t7_unmapped

    def _combine_alignments(self, sb_alignments, t7_alignments, bc_alignments):
        sel_columns = ['query_name', 'query_start', 'query_end', 'type']

        combined = sb_alignments[sel_columns].merge(t7_alignments[sel_columns],
                                                    on='query_name', how='inner', suffixes=['_SB', '_T7'])

        combined = combined.merge(bc_alignments[sel_columns + ['target_name', 'query_seq']],
                                  on='query_name', how='inner')

        return combined

    def _genomic_sequences(self, combined_alignment):
        query_names, query_seqs = combined_alignment['query_name'], combined_alignment['query_seq']
        sb_ends, t7_starts = combined_alignment['query_end_SB'], combined_alignment['query_start_T7']

        genomic_seqs = []
        for i, query_seq in enumerate(query_seqs):
            genomic_seq = query_seq[sb_ends.iloc[i]:t7_starts.iloc[i]]
            genomic_seqs.append(Sequence(query_names.iloc[i], genomic_seq))

        return genomic_seqs


class SBAlignmentStats(object):

    def __init__(self, reads, result):
        self.reads = reads
        self.result = result

    def plot_mapped(self):
        num_mapped = len(self.result.genomic_sequences)
        counts = pandas.Series({'mapped': num_mapped, 'unmapped': len(self.reads) - num_mapped })
        return ggplot_bar_hist(counts, xlabel='Mapping', title='Vector mapping (all)')

    def plot_unmapped_type(self, kind='pie'):
        sb_ids = set([s.name for s in self.result.vector_unmapped['SB']])
        t7_ids = set([s.name for s in self.result.vector_unmapped['T7']])
        both_ids = sb_ids.intersection(t7_ids)

        counts = pandas.Series({ 'No SB': len(self.result.vector_unmapped['SB']),
                                 'No T7': len(self.result.vector_unmapped['T7']),
                                 'No BC': len(self.result.barcode_unmapped),
                                 'No T7 and SB': len(both_ids) } )
        if kind == 'pie':
            plot = ggplot_pie(counts)
        elif kind == 'bar':
            plot = ggplot_bar_hist(counts, title='Unmapped cause', xlabel='Cause')
        else: raise ValueError

        return plot

    def plot_genomic_lengths(self):
        seq_lengths = [len(s.seq) for s in self.result.genomic_sequences]
        len_frame = pandas.DataFrame(seq_lengths, columns='x')
        return ggplot_hist(len_frame, binwidth=10)


    def plot_alignment_types(self):
        def type_count_frame(alignments, name):
            type_counts = alignments['type'].value_counts()
            count_frame = type_counts.div(float(type_counts.sum())).reset_index()
            count_frame.columns = ['type', 'count']
            count_frame['name'] = name
            return count_frame

        t7_counts = type_count_frame(self.result.vector_alignments['T7'], 'T7')
        sb_counts = type_count_frame(self.result.vector_alignments['SB'], 'SB')
        #bc_counts = type_count_frame(self.result.barcode_alignments, 'BC')
        count_frame = pandas.concat([t7_counts, sb_counts], ignore_index=True)

        return ggplot_stacked_bar(count_frame, xvar='name', yvar='count', fillvar='type',
                                  title='Alignment type', xlabel='Vector', ylabel='Fraction')

def _alignments_to_frame(alns):
    if len(alns) == 0: return None
    values = alns.values() if type(alns) == dict else alns
    return pandas.DataFrame(values, columns=values[0]._fields)