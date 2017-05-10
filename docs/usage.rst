=====
Usage
=====

Identifying insertions
----------------------

The **pyim-align** command is used to identify insertions using sequence reads
from targeted DNA-sequencing of insertion sites. The command provides access
to various pipelines which (in essence) perform the following functions:

    - Reads are filtered to remove reads that do not contain the correct
      technical sequences (such as transposon sequences or required adapter
      sequences).
    - Reads are trimmed to remove any non-genomic sequences (including
      transposon/adapter sequences and any other technical sequences). Reads
      that are too short after trimming are removed from the analysis, to
      avoid issues during alignment.
    - The remaining (genomic) reads are aligned to the reference genome.
    - The resulting alignment is analyzed to identify the location and
      orientation of the corresponding insertion sites.

The exact implementation of these steps differs between pipelines and depends
on the design of the sequencing experiment. For an overview of the available
pipelines, see :ref:`api_pipelines`.

An example of calling ``pyim-align`` using the ``shearsplink`` pipeline is
as follows:

.. code-block:: bash

    pyim-align shearsplink
        --reads ./reads.fastq.gz
        --bowtie_index /path/to/index
        --output_dir ./out
        --transposon /path/to/transposon.fa
        --linker /path/to/linker.fa

This produces an ``insertions.txt`` file in the ``./out`` directory,
describing the identified insertions.

Merging/splitting datasets
--------------------------

The **pyim-merge** command can be used to merge different sets of insertions.
This is mainly useful for combining insertions from multiple samples or from
different sequencing datasets. The basic command is as follows:

.. code-block:: bash

    pyim-merge --insertions ./sample1/insertions.txt \
                            ./sample2/insertions.txt \
               --output ./merged.txt

This command adds an additional ``sample`` column to the merged insertions
if this column is not yet present in the source files. By default, the used
sample names are derived from the names of the folders containing the
source files (``sample1`` and ``sample2`` in this example). These names can be
overridden using the ``--sample_names`` parameter.

Alternatively, the **pyim-split** command can be used to split a merged
insertion file (containing multiple samples) to obtain separate insertion
files for each sample. The basic command is as follows:

.. code-block:: bash

    pyim-split --insertions ./merged.txt \
               --output_dir ./out

A specific subset of samples can be extracted using the ``--samples`` argument.

Annotating insertions
---------------------

.. code-block:: bash

    pyim-annotate window --insertions ./out/insertions.txt
                         --output ./out/insertions.ann.txt
                         --gtf reference.gtf
                         --window_size 20000

Identifying CISs
----------------

TODO
