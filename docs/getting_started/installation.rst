.. _Miniconda: https://conda.pydata.org/miniconda.html
.. _Miniforge: https://github.com/conda-forge/miniforge
.. _Mamba: https://prefix.dev/docs/mamba
.. _Conda: https://conda.pydata.org
.. _PyPi: https://pypi.org/project/snakemake/
.. _Pixi: https://pixi.sh
.. _Micromamba: https://prefix.dev/docs/mamba
.. _WSL: https://learn.microsoft.com/en-us/windows/wsl/install


.. _getting_started-installation:

============
Installation
============

Snakemake is available through :ref:`pixi<pixi-install>`, :ref:`conda<conda-install>`, or :ref:`pip<pip-install>`.
Alternatively, you can also build it from source.

.. _pixi-install:

These instructions set out how to obtain and install the software and data on Linux. It is assumed that you have:

- access to the Bash or Zsh shell on a fairly modern Linux or macOS system
- sufficient disk space (~1GB) to store the software

You do **not** need root/administrator access.

.. note::
    Snakemake is intended to be run on Linux operating systems.
    If you want to run Snakemake on Windows, we recommend using the Windows Subsystem for Linux (WSL_), and then following the instructions for Linux.


Install via Pixi
-----------------

Pixi_ is a fast package management tool that acts as a drop-in replacement for Conda.

This is the **recommended** way to install Snakemake,
and will also enable Snakemake to use conda to :ref:`isolate your
workflow environments<integrated_package_management>`.

To install `pixi`, you can run the following command in your terminal:

.. code:: console

    $ curl -fsSL https://pixi.sh/install.sh | bash


If you wish to install snakemake globally on your system, run:

.. code-block:: console

    $ pixi global install snakemake -c conda-forge -c bioconda


You can also keep the installation contained in a pixi environment.
In the directory where you want to install snakemake, run:

.. code-block:: console

    $ pixi init; pixi add python; pixi add snakemake --pypi

You can then enter the environment with

.. code-block:: console

    $ pixi shell


.. _conda-install:

Install via Conda/Mamba
-----------------------

Conda_ is a language-agnostic package and environment manager.

There are many different ways to install conda; We recommend Miniforge_.

You can use the command below to download and install Miniconda, on your operating system:

.. tab-set::

    .. tab-item:: Linux

        .. code:: console

            $ curl -L https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh -o Miniforge3-Linux-x86_64.sh
            $ bash Miniforge3-Linux-x86_64.sh

    .. tab-item:: MacOS (x86_64)

        .. code:: console

            $ curl -L https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-MacOSX-x86_64.sh -o Miniforge3-MacOSX-x86_64.sh
            $ bash Miniforge3-MacOSX-x86_64.sh

    .. tab-item:: MacOS (arm64)

        .. code:: console

            $ curl -L https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-MacOSX-arm64.sh -o Miniforge3-MacOSX-arm64.sh
            $ bash Miniforge3-MacOSX-arm64.sh


Of course, any other conda-based package manager can be used as well, e.g. Pixi_, Mamba_ or Micromamba_.
Note however that for the :ref:`conda integration <integrated_package_management>` of Snakemake, for now it requires the `conda` command to be installed in the root environment or in the same environment as Snakemake itself.

Full Installation
-----------------

Snakemake can be installed with all the goodies needed to run in any environment and for creating interactive reports via

.. code-block:: console

    $ conda create -c conda-forge -c bioconda -c nodefaults -n snakemake snakemake

from the `Bioconda <https://bioconda.github.io>`_ channel.
This will install snakemake into an isolated software environment, that has to be activated with

.. code-block:: console

    $ conda activate snakemake
    $ snakemake --help

Installing into isolated environments is best practice in order to avoid side effects with other packages.

Minimal Installation
--------------------

A minimal version with only the necessary requirements can be installed with

.. code-block:: console

    $ conda create -c conda-forge -c bioconda -c nodefaults -n snakemake snakemake-minimal

Notes on Bioconda as a Package Source
-------------------------------------

Note that Snakemake is available via Bioconda for historical, reproducibility, and continuity reasons (although it is not limited to biology applications at all).
However, it is easy to combine Snakemake installation with other channels, e.g., by prefixing the package name with ``::bioconda``, i.e.,

.. code-block:: console

    $ conda activate base
    $ conda create -n some-env -c conda-forge -c nodefaults bioconda::snakemake ...


.. _pip-install:

Install via Pip
====================

pip is the default package manager in python.

You can install snakemake via pip with:

.. code-block:: console

    $ pip install snakemake


Install a Development Version
=============================================

If you want to quickly try out an unreleased version from the snakemake repository, e.g. to check whether a bug fix works for your workflow, you can get the current state of the main branch with:

.. code-block:: console

    $ conda create --only-deps -n snakemake-dev -c conda-forge -c bioconda -c nodefaults snakemake
    $ conda activate snakemake-dev
    $ pip install git+https://github.com/snakemake/snakemake

You can also install the current state of another branch or the repository state at a particular commit.
For information on the syntax for this, see `the pip documentation on git support <https://pip.pypa.io/en/stable/topics/vcs-support/#git>`_.


Editor Integrations
===================

There are a number of integrations with popular editors that offer convenient features such as syntax highlighting and formatting:

* `VSCode <https://github.com/snakemake/snakemake-lang-vscode-plugin>`_
* `Vim <https://github.com/snakemake/snakemake/tree/main/misc/vim>`_
* `Zed <https://github.com/lvignoli/zed-snakemake>`_
