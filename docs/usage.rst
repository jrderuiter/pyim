=====
Usage
=====

Identifying insertions
----------------------

Overview
========

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
on the design of the sequencing experiment.

Each pipeline takes...

Pipelines
=========

ShearSplink
~~~~~~~~~~~

The ``shearsplink`` pipeline is designed to analyze data from samples that
were sequenced using the ShearSplink_ protocol. ShearSplink sequence reads are
expected to have the following structure::

    [Transposon][Genomic][Linker]

Here, the ``transposon`` element represents part of the transposon sequence,
which is used (a) verify that the read does indeed involve the transposon and
(b) to determine the exact breakpoint between the transposon sequence and
the flanking genomic sequence (the ``genomic`` element in the sequence). The
``linker`` represents an adapter that is ligated to the (sheared) genomic as
part of the protocol. The position of the linker sequence allows us to assess
the depth/clonality of individual insertions by determining the number of
unique ligation points between the adapter and genomic DNA (see the
ShearSplink_ publication for more details).

The pipeline can be run using the following basic command:

.. code-block:: bash

    pyim-align --reads ./reads/Pool3.1.TCA.454Reads.fna \
               --output_dir ./out \
               --transposon ~/Software/python/pyim/data/sb.transposon.fa \
               --linker ~/Software/python/pyim/data/sb.linker.fa \
               --contaminants ~/Software/python/pyim/data/sb.contaminants.fa \
               --bowtie_index ~/References/mus_musculus/mm10/indices/bowtie2/Mus_musculus.GRCm38.dna.primary_assembly

The ``--linker`` and ``--contaminants`` arguments are optional. If the linker
sequence is omitted, the pipeline assumes that reads were sequenced without
including the linker. If a contaminants file is provided, the sequence reads
are first filtered for the provided contaminant sequences before further
processing. This enables filtering for specific contaminants, such as reads
stemming from the transposon donor locus.

.. _ShearSplink: https://www.ncbi.nlm.nih.gov/pubmed/21852388

Multiplexed ShearSplink
~~~~~~~~~~~~~~~~~~~~~~~

The ``shearsplink-multiplexed`` pipeline is an extended version of the
ShearSplink pipeline, which can handle multiplexed datasets. The main advantage
of the pipeline is that it directly tags insertions belonging to specific
samples, rather than us having to first demultiplex the sequence reads and then
having to analyze each sample individually.

The pipeline takes the same arguments as the basic ShearSplink pipeline, but
adds two arguments ``--barcodes`` and ``--barcode_mapping``. These arguments
specify which barcode sequences have been used to index reads and allow us
to specify an (optional) mapping from barcodes to sample names. If no mapping
is provided, insertions are tagged with the name of the corresponding barcode.

Multiplexed ShearSplink reads are expected to have the following structure::

    [Barcode][Transposon][Genomic][Linker]

The ``transposon``, ``genomic`` and ``linker`` elements are the same as for
the normal ShearSplink pipeline. The ``barcode`` sequence indicates from which
sample the read originated. Barcode sequences should correspond with a sequence
in the provided barcode file.

Nextera
~~~~~~~

TODO

Merging/splitting datasets
--------------------------

.. code-block:: bash

    pyim-merge --insertions ./out1/insertions.txt ./out2/insertions.txt \
               --output ./merged.txt

Annotating insertions
---------------------

.. code-block:: bash

    pyim-annotate window --insertions ./out/insertions.txt
                         --output ./out/insertions.ann.txt
                         --gtf reference.gtf
                         --window_size 20000

Identifying CISs
----------------

