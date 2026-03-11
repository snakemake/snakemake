
.. _tutorial-setup:

Setup
====

.. _Snakemake: https://snakemake.readthedocs.io
.. _Snakemake homepage: https://snakemake.readthedocs.io
.. _GNU Make: https://www.gnu.org/software/make
.. _Python: https://www.python.org
.. _BWA: http://bio-bwa.sourceforge.net
.. _SAMtools: https://www.htslib.org
.. _BCFtools: https://www.htslib.org
.. _Pandas: https://pandas.pydata.org
.. _Miniconda: https://conda.pydata.org/miniconda.html
.. _Miniforge: https://github.com/conda-forge/miniforge
.. _Mamba: https://github.com/mamba-org/mamba
.. _Conda: https://conda.pydata.org
.. _Pixi: https://pixi.sh/
.. _Pixi installation: https://pixi.sh/latest/#installation
.. _Pixi automated switching: https://pixi.sh/latest/switching_from/conda/#automated-switching
.. _Bash: https://www.tldp.org/LDP/Bash-Beginners-Guide/html
.. _Atom: https://atom.io
.. _Graphviz: https://www.graphviz.org
.. _PyYAML: https://pyyaml.org
.. _Docutils: https://docutils.sourceforge.io
.. _Jinja2: https://jinja.palletsprojects.com
.. _NetworkX: https://networkx.github.io
.. _Matplotlib: https://matplotlib.org
.. _Pysam: https://pysam.readthedocs.io
.. _Bioconda: https://bioconda.github.io
.. _WSL: https://docs.microsoft.com/en-us/windows/wsl/about
.. _WSL Documentation: https://docs.microsoft.com/en-us/windows/wsl/install-win10
.. _Vagrant: https://www.vagrantup.com
.. _Vagrant Documentation: https://docs.vagrantup.com
.. _Blogpost: https://blog.osteel.me/posts/2015/01/25/how-to-use-vagrant-on-windows.html


.. _tutorial-free-on-gitpod:

Run tutorial for free in the cloud via Gitpod
----

The easiest way to run this tutorial is to use Gitpod, which enables performing the exercises via your browser---including all required software, for free and in the cloud.
In order to do this, simply open the predefined `snakemake-tutorial GitPod workspace <https://gitpod.io/#https://github.com/snakemake/snakemake-tutorial-data>`_ in your browser.
GitPod provides you with a `Theia development environment <https://theia-ide.org/docs>`_, which you can learn about in the linked documentation.
Once you have a basic understanding of this environment, you can go on directly with :ref:`tutorial-basics`.

.. note::
    A common thing to happen while using the development environment in GitPod is to hit ``Ctrl-s`` while in the terminal window, because you wanted to save a file in the editor window.
    This will freeze your terminal.
    To get it back, make sure you selected the terminal window by clicking on it and then hit ``Ctrl-q``.

Running the tutorial on your local machine
----

If you prefer to run the tutorial on your local machine, please follow the installation instructions
for either pixi or miniforge.

The tutorial assumes that you are using either Linux or MacOS X.
Snakemake, Pixi_, and Miniforge_ work also under Windows, but the Windows shell is too different to be able to provide generic examples.


1. Install Snakemake
::::::::::::::::::::

Follow the :ref:`installation instructions<getting_started-installation>` using either pixi or miniforge options to install snakemake.


2. Prepare a working directory
:::::::::::::::::::::::::::::::::::::

Create a new directory ``snakemake-tutorial`` at a **place you can easily remember** and change into that directory in your terminal:

.. code:: console

    $ mkdir snakemake-tutorial
    $ cd snakemake-tutorial

In this directory, we will later create an example workflow that illustrates the Snakemake syntax and execution environment.

Next, we will download some example data on which to run the workflow:

.. code:: console

    $ curl -L https://api.github.com/repos/snakemake/snakemake-tutorial-data/tarball -o snakemake-tutorial-data.tar.gz

Extract the data:

.. tab-set::

    .. tab-item:: Linux

        .. code:: console

            $ tar --wildcards -xf snakemake-tutorial-data.tar.gz --strip 1 "*/data" "*/environment.yaml"

    .. tab-item:: MacOS (x86_64)

        .. code:: console

            $ tar -xf snakemake-tutorial-data.tar.gz --strip 1 "*/data" "*/environment.yaml"

This will create a folder ``data`` and a file ``environment.yaml`` in the working directory.

3. Install Workflow Dependencies
::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

The tutorial provides an ``environment.yaml`` file that lists all required software for the workflow.

To import this environment and install all the required software for the tutorial workflow:

.. tab-set::
    :sync-group: package-manager

    .. tab-item:: Pixi
        :sync: pixi

        .. code:: console

            $ pixi init --import environment.yaml

    .. tab-item:: Conda
        :sync: conda

        .. code:: console

            $ conda activate base
            $ mamba env create --name snakemake-tutorial --file environment.yaml

4. Activate the environment
::::::::::::::::::::::::::::::::::

To activate the ``snakemake-tutorial`` environment, execute

.. tab-set::
    :sync-group: package-manager

    .. tab-item:: Pixi
        :sync: pixi

        .. code:: console

            $ pixi shell

    .. tab-item:: Conda
        :sync: conda

        .. code:: console

            $ conda activate snakemake-tutorial

Now you can use the installed tools within this shell.

If you want to exit this shell later, use

.. tab-set::
    :sync-group: package-manager

    .. tab-item:: Pixi
        :sync: pixi

        .. code:: console

            $ exit

    .. tab-item:: Conda
        :sync: conda

        .. code:: console

            $ conda deactivate

to return to your normal shell.
But **don't do that now**, since we finally want to start working with Snakemake! :-)
