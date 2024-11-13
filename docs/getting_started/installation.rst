.. _Miniconda: https://conda.pydata.org/miniconda.html
.. _Miniforge: https://github.com/conda-forge/miniforge
.. _Mamba: https://github.com/mamba-org/mamba
.. _Conda: https://conda.pydata.org
.. _PyPi: https://pypi.org/project/snakemake/


.. _getting_started-installation:

============
Installation
============

Snakemake is available on PyPi_ as well as through Conda_ and also from source code.
You can use one of the following ways for installing Snakemake.

.. _conda-install:

Installation via Conda/Mamba
============================

This is the **recommended** way to install Snakemake,
because it also enables Snakemake to :ref:`handle software dependencies of your
workflow <integrated_package_management>`.

First, you need to install a Conda-based Python3 distribution.
The recommended choice is Miniforge_.

Of course, any other conda-based package manager can be used as well, e.g. Pixi_, Mamba_ or Micromamba_.
Note however that for the :ref:`conda integration <integrated_package_management>` of Snakemake, for now it requires the `conda` command to be installed in the root environment or in the same environment as Snakemake itself.

Full installation
-----------------

Snakemake can be installed with all goodies needed to run in any environment and for creating interactive reports via

.. code-block:: console

    $ conda create -c conda-forge -c bioconda -n snakemake snakemake

from the `Bioconda <https://bioconda.github.io>`_ channel.
This will install snakemake into an isolated software environment, that has to be activated with

.. code-block:: console

    $ conda activate snakemake
    $ snakemake --help

Installing into isolated environments is best practice in order to avoid side effects with other packages.

Minimal installation
--------------------

A minimal version with only the necessary requirements can be installed with

.. code-block:: console

    $ conda create -c conda-forge -c bioconda -n snakemake snakemake-minimal

Notes on Bioconda as a package source
-------------------------------------

Note that Snakemake is available via Bioconda for historical, reproducibility, and continuity reasons (although it is not limited to biology applications at all).
However, it is easy to combine Snakemake installation with other channels, e.g., by prefixing the package name with ``::bioconda``, i.e.,

.. code-block:: console

    $ conda activate base
    $ conda create -n some-env -c conda-forge bioconda::snakemake ...

Installation via pip
====================

Instead of conda, snakemake can be installed with pip:

.. code-block:: console

    $ pip install snakemake


Installation of a development version via pip
=============================================

If you want to quickly try out an unreleased version from the snakemake repository (which you cannot get via e.g. bioconda or PyPi_, yet), for example to check whether a bug fix works for you workflow, you can get the current state of the main branch with:

.. code-block:: console

    $ conda create --only-deps -n snakemake-dev snakemake
    $ conda activate snakemake-dev
    $ pip install git+https://github.com/snakemake/snakemake

You can also install the current state of another branch or the repository state at a particular commit.
For information on the syntax for this, see `the pip documentation on git support <https://pip.pypa.io/en/stable/topics/vcs-support/#git>`_.


Editor integrations
===================

* `VSCode <https://github.com/snakemake/snakemake-lang-vscode-plugin>`_
* `Vim <https://github.com/snakemake/snakemake/tree/main/misc/vim>`_
* `Zed <https://github.com/lvignoli/zed-snakemake>`_
