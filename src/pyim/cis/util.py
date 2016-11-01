import itertools
import operator

from frozendict import frozendict
import numpy as np
import toolz


def assign_strand(cis_sites, insertions, mapping, min_homogeneity=0.75):
    """Assigns CIS sites the average strand of their insertions."""

    ins_lookup = {insertion.id: insertion for insertion in insertions}

    for cis_site in cis_sites:
        # Lookup strands of CIS insertions.
        cis_strands = np.array([ins_lookup[ins_id].strand
                                for ins_id in mapping[cis_site.id]])

        # Calculate mean strand, strand and homogeneity.
        mean_strand = np.mean(cis_strands)
        strand = np.sign(mean_strand)
        homogeneity = np.sum(cis_strands == strand) / len(cis_strands)

        # If homogeneity is below the given threshold, then we don't
        # assign a specific strand (signified by a nan).
        if homogeneity < min_homogeneity:
            strand = np.nan

        # Merge strand metadata with existing metadata.
        strand_metadata = {'strand_mean': mean_strand,
                           'strand_homogeneity': homogeneity}
        metadata = toolz.merge(cis_site.metadata, strand_metadata)

        yield cis_site._replace(strand=strand, metadata=frozendict(metadata))


def invert_otm_mapping(mapping):
    """Inverts a one-to-many mapping."""

    # Create a list of inverted (v, k) tuples.
    tuples = (itertools.zip_longest(v, [k], fillvalue=k)
              for k, v in mapping.items() if len(v) > 0) # yapf: disable
    tuples = itertools.chain.from_iterable(tuples)

    # Sort tuples by first element.
    id_attr = operator.itemgetter(0)
    sorted_tuples = sorted(tuples, key=id_attr)

    # Create inverted dictionary using groupby.
    inverted = {k: set(list(zip(*grp))[1])
                for k, grp in itertools.groupby(
                    sorted_tuples, key=id_attr)}

    return inverted
