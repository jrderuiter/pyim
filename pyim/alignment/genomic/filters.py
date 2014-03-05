

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
