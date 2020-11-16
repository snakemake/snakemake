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

The default conda solver is a bit slow and sometimes has issues with `selecting the latest package releases <https://github.com/conda/conda/issues/9905>`_. Therefore, we recommend to install `Mamba <https://github.com/mamba-org/mamba>`_ as a drop-in replacement via

.. code-block:: console

    $ conda install -c conda-forge mamba

Full installation
-----------------

Snakemake can be installed with all goodies needed to run in any environment and for creating interactive reports via

.. code-block:: console

    $ mamba create -c conda-forge -c bioconda -n snakemake snakemake

from the `Bioconda <https://bioconda.github.io>`_ channel.
This will install snakemake into an isolated software environment, that has to be activated with

.. code-block:: console

    $ conda activate snakemake
    $ snakemake --help

Installing into isolated environments is best practice in order to avoid side effects with other packages.

Note that full installation is not possible from **Windows**, because some of the dependencies are Unix (Linux/MacOS) only.
For Windows, please use the minimal installation below.

Minimal installation
--------------------

A minimal version of Snakemake which only depends on the bare necessities can be installed with

.. code-block:: console

    $ mamba create -c bioconda -c conda-forge -n snakemake snakemake-minimal

In contrast to the full installation, which depends on some Unix (Linux/MacOS) only packages, this also works on Windows.

Notes on Bioconda as a package source
-------------------------------------

Note that Snakemake is available via Bioconda for historical, reproducibility, and continuity reasons (although it is not limited to biology applications at all).
However, it is easy to combine Snakemake installation with other channels, e.g., by prefixing the package name with ``::bioconda``, i.e.,

.. code-block:: console

    $ mamba create -n some-env -c conda-forge bioconda::snakemake bioconda::snakemake-minimal ...

Installation via pip
====================

Instead of conda, snakemake can be installed with pip.
However, note that snakemake has non-python dependencies, such that the pip based installation has a limited functionality if those dependencies are not manually installed in addition.
