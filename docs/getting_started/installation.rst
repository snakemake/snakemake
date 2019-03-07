.. getting_started-installation:

============
Installation
============

Snakemake is available on PyPi as well as through Bioconda and also from source code.
You can use one of the following ways for installing Snakemake.

Installation via Conda
======================

This is the recommended way to install Snakemake,
because it also enables Snakemake to :ref:`handle software dependencies of your
workflow <integrated_package_management>`.

First, you have to install the Miniconda Python3 distribution.
See `here <https://conda.io/docs/install/quick.html>`_ for installation instructions.
Make sure to ...

* Install the **Python 3** version of Miniconda.
* Answer yes to the question whether conda shall be put into your PATH.

Then, you can install Snakemake with

.. code-block:: console

    $ conda install -c bioconda -c conda-forge snakemake

from the `Bioconda <https://bioconda.github.io>`_ channel.
A minimal version of Snakemake which only depends on the bare necessities can be installed with

.. code-block:: console

    $ conda install -c bioconda -c conda-forge snakemake-minimal

Note that Snakemake is available via Bioconda for historical, reproducibility, and continuity reasons.
However, it is easy to combine Snakemake installation with other channels, e.g., by prefixing the package name with ``::bioconda``, i.e.,

.. code-block:: console

    $ conda install -c conda-forge bioconda::snakemake bioconda::snakemake-minimal

Global Installation
===================

With a working Python ``>=3.5`` setup, installation of Snakemake can be performed by issuing

.. code-block:: console

    $ easy_install3 snakemake

or

.. code-block:: console

    $ pip3 install snakemake

in your terminal.


Installing in Virtualenv
========================

To create an installation in a virtual environment, use the following commands:

.. code-block:: console

    $ virtualenv -p python3 .venv
    $ source .venv/bin/activate
    $ pip install snakemake


Installing from Source
======================

We recommend installing Snakemake into a virtualenv instead of globally.
Use the following commands to create a virtualenv and install Snakemake.
Note that this will install the development version and as you are installing from the source code, we trust that you know what you are doing and how to checkout individual versions/tags.

.. code-block:: console

    $ git clone https://bitbucket.org/snakemake/snakemake.git
    $ cd snakemake
    $ virtualenv -p python3 .venv
    $ source .venv/bin/activate
    $ python setup.py install

You can also use ``python setup.py develop`` to create a "development installation" in which no files are copied but a link is created and changes in the source code are immediately visible in your ``snakemake`` commands.
