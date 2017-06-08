.. highlight:: shell

============
Installation
============

Dependencies
------------

PyIM requires Python 3.4 and has been tested on macOS and Linux.

The following external dependencies are also required for full functionality:

- Bowtie2
- Cutadapt
- CIMPL (R package, via rpy2)

These external tools should be available in ``$PATH``. CIMPL, which is an R
package, should be loadable in the default R installation.

Using bioconda (recommended)
----------------------------

The recommended way to install PyIM is using bioconda, as this installs
PyIM together with its external dependencies into an isolated environment
using a single command:

.. code:: bash

  conda create -n pyim pyim

This assumes that condas channels have been configured as recommended_
by bioconda.

.. _recommended: https://bioconda.github.io/#set-up-channels

Alternatively, PyIM can be installed in an existing environment using:

.. code:: bash

  conda install pyim

Using pip
---------

IM-Fusion can also be installed from PyPI using pip:

.. code:: bash

    pip install pyim

Note that this does not install any of the required external dependencies,
which must therefore be installed manually.
