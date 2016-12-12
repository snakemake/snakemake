.. getting_started-installation:

============
Installation
============

Snakemake is available on PyPi as well as through Bioconda and also from source code.
You can use one of the following ways for installing Snakemake.


Global Installation
===================

With a working Python 3 setup, installation of Snakemake can be performed by issuing

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


Installing Conda
================

In case you have to install Python 3 yourself, we recommend to use the Miniconda Python 3 distribution (http://conda.pydata.org/miniconda.html).

With Miniconda installed, you can issue

.. code-block:: console

    $ conda install -c bioconda snakemake

to install Snakemake from the bioconda channel.


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
