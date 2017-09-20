from pyim.align.util import AlignmentSummary


class TestAlignmentSummary(object):
    """Unit tests for the AlignmentSummary class."""

    def test_merge_within_distance(self):
        """Test basic call to merge_within_distance.

        Should merge entries around 58654912-58654918.
        """

        summary = AlignmentSummary(values={
            's1': {
                ('6', 30484510, 1): [30484517],
                ('6', 56984262, 1): [56984272],
                ('6', 58654912, 1): [58654922],
                ('6', 58654914, 1): [58654922],
                ('6', 58654916, 1): [58654923],
                ('6', 58654918, 1): [58654925],
                ('6', 102007479, 1): [102007489],
                ('6', 139746596, 1): [139746606],
                ('6', 58654916, -1): [58654926]
            }
        })

        # Merge within distance of 10.
        summary_merged = summary.merge_within_distance(max_dist=10)

        # Check merge result.
        assert summary_merged.values == {
            's1': {
                ('6', 58654916, -1): [58654926],
                ('6', 30484510, 1): [30484517],
                ('6', 102007479, 1): [102007489],
                ('6', 139746596, 1): [139746606],
                ('6', 58654915, 1): [58654922, 58654922, 58654923, 58654925],
                ('6', 56984262, 1): [56984272]
            }
        }

    def test_merge_within_distance_zero(self):
        """Test call to merge_within_distance with max_dist = 0.

        Should leave values unchanged.
        """

        values = {
            's1': {
                ('6', 30484510, 1): [30484517],
                ('6', 56984262, 1): [56984272],
                ('6', 58654912, 1): [58654922],
                ('6', 58654914, 1): [58654922],
                ('6', 58654916, 1): [58654923],
                ('6', 58654918, 1): [58654925],
                ('6', 102007479, 1): [102007489],
                ('6', 139746596, 1): [139746606],
                ('6', 58654916, -1): [58654926]
            }
        }

        summary = AlignmentSummary(values=values)
        summary_merged = summary.merge_within_distance(max_dist=0)

        assert summary_merged.values == values
