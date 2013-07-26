
import MySQLdb
import pandas.io.sql

from collections import defaultdict
from bx.intervals.intersection import IntervalTree, Interval


class EnsemblGeneDB(object):

    def __init__(self):
        db = MySQLdb.connect(host="ensembldb.ensembl.org", port=5306,
                             user="anonymous", db="mus_musculus_core_67_37")
        self._frame = pandas.io.sql.read_frame('SELECT * FROM gene', db)

    def _create_interval_tree(self):
        region = self._frame['seq_region_id'].values
        start = self._frame['seq_region_start'].values
        end = self._frame['seq_region_end'].values
        strand = self._frame['seq_region_strand'].values

        trees = defaultdict(IntervalTree)
        for i in range(len(self._frame)):
            interval = Interval()
            trees[region[i]].add(interval)
