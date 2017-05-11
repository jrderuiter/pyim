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

To install PyIM, run this command in your terminal:

.. code-block:: console

    $ pip install https://github.com/jrderuiter/pyim/archive/0.2.0.tar.gz

If you don't have `pip`_ installed, this `Python installation guide`_ can guide
you through the process.

.. _pip: https://pip.pypa.io
.. _Python installation guide: http://docs.python-guide.org/en/latest/starting/installation/

License
-------

This software is released under the MIT license.
