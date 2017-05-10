===============================
PyIM
===============================

.. image:: https://img.shields.io/travis/jrderuiter/pyim.svg
        :target: https://travis-ci.org/jrderuiter/pyim

PyIM (Python Insertional Mutagenesis) is a python package for analyzing
insertional mutagenesis data from targeted sequencing of transposon insertion
sites. The package provides several command line tools for identifying
insertions, calling common insertion sites (CISs) and annotating
insertions/CISs directly from the command line. It also aims to provides
the basic building blocks for implementing new pipelines, CIS callers, etc.

Documentation
-------------

PyIM's documentation is available at
`jrderuiter.github.io/pyim <http://jrderuiter.github.io/pyim/>`_.


Requirements
------------

PyIM is written for Python 3 and requires Python 3.3 or newer to be installed.
Depending on the used functionality, PyIM also has the following external
dependencies:

- cutadapt/bowtie2 (Needed for identifying insertions from sequencing data)
- cimpl (R package, needed for calling CIS sites using CIMPL)

Installation
------------

Using conda
~~~~~~~~~~~

The recommended way to install PyIM is using conda, as with conda you can
install PyIM together with its external dependencies (cutadapt and bowtie2)
into an isolated environment using a single command:

.. code:: bash

    conda create -n pyim -c jrderuiter -c bioconda -c r pyim

Alternatively, PyIM can be installed in an existing environent using:

.. code:: bash

    conda install -c jrderuiter -c bioconda -c r pyim

Conda packages are available for both OSX and Linux (64-bit).

Using pip
~~~~~~~~~

PyIM can be installed from Github using pip as follows:

.. code:: bash

    pip install git+git://github.com/jrderuiter/pyim.git#egg=pyim

Note that in this case, external dependencies must be installed manually.

Unfortunately, PyIM is not yet available on PyPI, though this may
change when the package is further developed.

License
-------

This software is released under the MIT license.
