.. _Miniconda: https://conda.pydata.org/miniconda.html
.. _Mambaforge: https://github.com/conda-forge/miniforge#mambaforge
.. _Micromamba: https://mamba.readthedocs.io/en/latest/user_guide/micromamba.html
.. _Mamba: https://github.com/mamba-org/mamba
.. _Conda: https://conda.pydata.org


.. _getting_started-installation:

============
Installation
============

Snakemake is available on PyPi as well as through Bioconda and also from source code.
You can use one of the following ways for installing Snakemake.

.. _conda-install:

Installation via Conda/Mamba
============================

This is the **recommended** way to install Snakemake,
because it also enables Snakemake to :ref:`handle software dependencies of your
workflow <integrated_package_management>`.

First, you need to install a Conda-based Python3 distribution.
The recommended choice is Mambaforge_, which not only provides the required Python and Conda commands, 
but also includes Mamba_, an extremely fast and robust replacement for the Conda_ package manager which is highly recommended.
The default conda solver is a bit slow and sometimes has issues with `selecting the latest package releases <https://github.com/conda/conda/issues/9905>`_. 
Therefore, we recommend to in any case use Mamba_; in case you are looking for something a bit more lightweight, you can also try Micromamba_ (and just replace all `mamba` commands with `micromamba`).

In case you don't use Mambaforge_ you can always install Mamba_ or Micromamba_ into any other Conda-based Python distribution with

.. code-block:: console

    $ conda install -n base -c conda-forge mamba micromamba

Full installation
-----------------

Snakemake can be installed with all goodies needed to run in any environment and for creating interactive reports via

.. code-block:: console

    $ mamba create -c conda-forge -c bioconda -n snakemake snakemake

from the `Bioconda <https://bioconda.github.io>`_ channel.
This will install snakemake into an isolated software environment, that has to be activated with

.. code-block:: console

    $ mamba activate snakemake
    $ snakemake --help

Installing into isolated environments is best practice in order to avoid side effects with other packages.

Notes on Bioconda as a package source
-------------------------------------

Note that Snakemake is available via Bioconda for historical, reproducibility, and continuity reasons (although it is not limited to biology applications at all).
However, it is easy to combine Snakemake installation with other channels, e.g., by prefixing the package name with ``::bioconda``, i.e.,

.. code-block:: console

    $ mamba activate base
    $ mamba create -n some-env -c conda-forge bioconda::snakemake ...

Installation via pip
====================

Instead of conda, snakemake can be installed with pip.
However, note that snakemake has non-python dependencies, such that the pip based installation has a limited functionality if those dependencies are not manually installed in addition.

A list of Snakemake's dependencies can be found within its `meta.yaml conda recipe <https://bioconda.github.io/recipes/snakemake/README.html>`_.


Installation of a development version via pip
=============================================

If you want to quickly try out an unreleased version from the snakemake repository (which you cannot get via bioconda, yet), for example to check whether a bug fix works for you workflow, you can get the current state of the main branch with:

.. code-block:: console

    $ mamba create --only-deps -n snakemake-dev snakemake
    $ mamba activate snakemake-dev
    $ pip install git+https://github.com/snakemake/snakemake

You can also install the current state of another branch or the repository state at a particular commit.
For information on the syntax for this, see `the pip documentation on git support <https://pip.pypa.io/en/stable/topics/vcs-support/#git>`_.
