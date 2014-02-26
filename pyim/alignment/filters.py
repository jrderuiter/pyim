
import logging, pandas
from pyim.io import Sequence


class AlignmentFilter(object):

    def __init__(self, verbose=False):
        self.verbose = verbose
        self.logger = logging.getLogger(self.__class__.__name__)

    def mask(self, alignment):
        raise NotImplementedError

    def apply(self, alignment):
        mask = self.mask(alignment)

        passed = alignment[mask]
        failed = alignment[~mask]
        dropped = self._dropped(passed, failed)

        if self.verbose:
            self._print_stats(alignment, failed, dropped)

        return passed, failed, dropped

    def _dropped(self, passed, failed):
        dropped_names = set(failed['query_name']).difference(passed['query_name'])
        dropped_frame = failed[failed['query_name'].isin(dropped_names)]
        seq_frame = dropped_frame.groupby('query_name')['query_seq'].first()
        return [Sequence(*sr) for sr in zip(seq_frame.index, seq_frame.values)]

    def _print_stats(self, alignment, failed, dropped):
        failed_fraction = len(failed)/float(len(alignment))
        self.logger.info('%d alignments failed filter (%02.2f%%)' % (len(failed), failed_fraction * 100))
        self.logger.info('%d reads were dropped as a result' % len(dropped))


#class UniqueAlignmentFilter(AlignmentFilter):
#
#    def apply(self, alignment):
#        result = super(UniqueAlignmentFilter, self).apply(alignment)
#        if result['query_name'].duplicated().sum() > 0:
#            print 'WARNING: duplicate alignments returned'
#        return result


class QueryEndFilter(AlignmentFilter):

    def mask(self, alignment):
        query_len = alignment['query_seq'].apply(len)
        return alignment['query_end'] == query_len


class BestScoreFilter(AlignmentFilter):

    def mask(self, alignment):
        groups = alignment.groupby('query_name')
        return groups['score'].transform(max) == alignment['score']


class MismatchFilter(AlignmentFilter):

    def __init__(self, n_max, include_start=False, include_end=False, verbose=False):
        super(MismatchFilter, self).__init__(verbose)
        self.n_max = n_max
        self.include_start = include_start
        self.include_end = include_end

    def mask(self, alignment):
        n_mismatched = alignment['alignment'].str.count('S')
        if self.include_start: n_mismatched += alignment['target_start']
        if self.include_end: n_mismatched += alignment['target_end']
        return n_mismatched <= self.n_max


class CombinedFilter(AlignmentFilter):

    def __init__(self, filters, verbose=False):
        super(CombinedFilter, self).__init__(verbose)
        self.filters = filters

    def mask(self, alignment):
        mask = pandas.Series([True] * len(alignment), index=alignment.index)
        for filt in self.filters:
            mask = filt.mask(alignment) & mask
        return mask