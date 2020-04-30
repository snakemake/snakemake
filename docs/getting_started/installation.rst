.. _getting_started-installation:

============
Installation
============

Snakemake is available on PyPi as well as through Bioconda and also from source code.
You can use one of the following ways for installing Snakemake.

.. _conda-install:

Installation via Conda
======================

This is the **recommended** way to install Snakemake,
because it also enables Snakemake to :ref:`handle software dependencies of your
workflow <integrated_package_management>`.

First, you have to install the Miniconda Python3 distribution.
See `here <https://conda.io/en/latest/miniconda.html>`_ for installation instructions.
Make sure to ...

* Install the **Python 3** version of Miniconda.
* Answer yes to the question whether conda shall be put into your PATH.

Then, you can install Snakemake with

.. code-block:: console

    $ conda create -c conda-forge -c bioconda -n snakemake snakemake=|version|

from the `Bioconda <https://bioconda.github.io>`_ channel.
This will install snakemake into an isolated software environment, that has to be activated with

.. code-block:: console

    $ conda activate snakemake
    $ snakemake --help

Installing into isolated environments is best practice in order to avoid side effects with other packages.
A minimal version of Snakemake which only depends on the bare necessities can be installed with

.. code-block:: console

    $ conda create -c bioconda -c conda-forge -n snakemake snakemake-minimal=|version|

Note that Snakemake is available via Bioconda for historical, reproducibility, and continuity reasons.
However, it is easy to combine Snakemake installation with other channels, e.g., by prefixing the package name with ``::bioconda``, i.e.,

.. parsed-literal::

    $ conda create -n some-env -c conda-forge bioconda::snakemake=|version| bioconda::snakemake-minimal=|version| ...

Installation via pip
====================

Instead of conda, snakemake can be installed with pip.
However, note that snakemake has non-python dependencies, such that the pip based installation has a limited functionality if those dependencies are not manually installed in addition.
