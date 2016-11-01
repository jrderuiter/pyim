API
============

Alignment
---------------


Vector
~~~~~~~~~~~~

.. autoclass:: pyim.alignment.vector.Alignment
  :members:
.. autofunction:: pyim.alignment.vector.align_exact(target, query, query_strand=1)
.. autofunction:: pyim.alignment.vector.align_ssw(target, query, query_strand=1)
.. autofunction:: pyim.alignment.vector.align_with_reverse(target, query, align_func, query_strand=1, **kwargs)
.. autofunction:: pyim.alignment.vector.align_multiple(target, queries, align_func, raise_error=False, **kwargs)
.. autofunction:: pyim.alignment.vector.align_chained(align_chained(target, query, align_funcs, **kwargs)
.. autofunction:: pyim.alignment.vector.compose
.. autofunction:: pyim.alignment.vector.filter_and(target, query, align_func, filters, **kwargs)
.. autofunction:: pyim.alignment.vector.filter_or(target, query, align_func, filters, **kwargs)
.. autofunction:: pyim.alignment.vector.filter_score(alignment, min_score)
.. autofunction:: pyim.alignment.vector.filter_coverage(alignment, min_coverage, min_identity)
.. autofunction:: pyim.alignment.vector.filter_end_match(alignment)

Genome
~~~~~~~~~~~~

.. autofunction:: pyim.alignment.bowtie2.align

Annotation
---------------

Annotators
~~~~~~~~~~~~

.. autofunction:: pyim.annotation.annotator.annotate_windows
.. autoclass:: pyim.annotation.annotator.Window
  :members:
.. autofunction:: pyim.annotation.annotator.annotate_rbm
.. autofunction:: pyim.annotation.annotator.annotate_rbm_cis

Metadata
~~~~~~~~~~~~

.. autofunction:: pyim.annotation.metadata.add_metadata
.. autofunction:: pyim.annotation.metadata.feature_distance
.. autofunction:: pyim.annotation.metadata.feature_orientation

Filtering
~~~~~~~~~~~~

.. autofunction:: pyim.annotation.filtering.filter_blacklist
.. autofunction:: pyim.annotation.filtering.select_closest
