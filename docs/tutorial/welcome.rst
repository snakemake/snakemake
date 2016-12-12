.. _tutorial-welcome:

==============
Tutorial Setup
==============

.. _Snakemake: http://snakemake.bitbucket.org
.. _Snakemake homepage: http://snakemake.bitbucket.org
.. _GNU Make: https://www.gnu.org/software/make
.. _Python: http://www.python.org
.. _BWA: http://bio-bwa.sourceforge.net
.. _SAMtools: http://www.htslib.org
.. _BCFtools: http://www.htslib.org
.. _Pandas: http://pandas.pydata.org
.. _Miniconda: http://conda.pydata.org/miniconda.html
.. _Conda: http://conda.pydata.org
.. _Bash: http://www.tldp.org/LDP/Bash-Beginners-Guide/html
.. _Atom: https://atom.io
.. _Anaconda: https://anaconda.org
.. _Graphviz: http://www.graphviz.org
.. _RestructuredText: http://docutils.sourceforge.net/rst.html
.. _data URI: https://developer.mozilla.org/en-US/docs/Web/HTTP/data_URIs
.. _JSON: http://json.org
.. _YAML: http://yaml.org
.. _DRMAA: http://www.drmaa.org
.. _rpy2: http://rpy.sourceforge.net
.. _R: https://www.r-project.org
.. _Rscript: https://stat.ethz.ch/R-manual/R-devel/library/utils/html/Rscript.html
.. _PyYAML: http://pyyaml.org
.. _Docutils: http://docutils.sourceforge.net
.. _Bioconda: https://bioconda.github.io
.. _Vagrant: https://www.vagrantup.com
.. _Vagrant Documentation: https://docs.vagrantup.com
.. _Blogpost: http://blog.osteel.me/posts/2015/01/25/how-to-use-vagrant-on-windows.html
.. _slides: http://slides.com/johanneskoester/deck-1

This tutorial introduces the text-based workflow system Snakemake_.
Snakemake follows the `GNU Make`_ paradigm: workflows are defined in terms of rules that define how to create output files from input files.
Dependencies between the rules are determined automatically, creating a DAG (directed acyclic graph) of jobs that can be automatically parallelized.

Snakemake sets itself apart from existing text-based workflow systems in the following way.
Hooking into the Python interpreter, Snakemake offers a definition language that is an extension of Python_ with syntax to define rules and workflow specific properties.
This allows to combine the flexibility of a plain scripting language with a pythonic workflow definition.
The Python language is known to be concise yet readable and can appear almost like pseudo-code.
The syntactic extensions provided by Snakemake maintain this property for the definition of the workflow.
Further, Snakemakes scheduling algorithm can be constrained by priorities, provided cores and customizable resources and it provides a generic support for distributed computing (e.g., cluster or batch systems).
Hence, a Snakemake workflow scales without modification from single core workstations and multi-core servers to cluster or batch systems.

While the examples presented here come from Bioinformatics, Snakemake is considered a general-purpose workflow management system for any discipline.

Also have a look at the corresponding slides_.

.. _tutorial-setup:

Requirements
::::::::::::

To go through this tutorial, you need the following software installed:

* Python_ ≥3.3
* Snakemake_ 3.9.0
* BWA_ 0.7.12
* SAMtools_ 1.3.1
* BCFtools_ 1.3.1
* Graphviz_ 2.38.0
* PyYAML_ 3.11
* Docutils_ 0.12

The easiest way to setup these prerequisites is to use the Miniconda_ Python 3 distribution.
The tutorial assumes that you are using either Linux or MacOS X.
Both Snakemake and Miniconda work also under Windows, but the Windows shell is too different to be able to provide generic examples.

Setup a Linux VM with Vagrant under Windows
:::::::::::::::::::::::::::::::::::::::::::

If you already use Linux or MacOS X, go on with **Step 1**.
If you use Windows, you can setup a Linux virtual machine (VM) with Vagrant_.
First, install Vagrant following the installation instructions in the `Vagrant Documentation`_.
Then, create a reasonable new directory you want to share with your Linux VM, e.g., create a folder ``vagrant-linux`` somewhere.
Open a command line prompt, and change into that directory.
Here, you create a 64-bit Ubuntu Linux environment with

.. code:: console

    > vagrant init hashicorp/precise64
    > vagrant up

If you decide to use a 32-bit image, you will need to download the 32-bit version of Miniconda in the next step.
The contents of the ``vagrant-linux`` folder will be shared with the virtual machine that is set up by vagrant.
You can log into the virtual machine via

.. code:: console

    > vagrant ssh

If this command tells you to install an SSH client, you can follow the instructions in this Blogpost_.
Now, you can follow the steps of our tutorial from within your Linux VM.


Step 1: Installing Miniconda 3
::::::::::::::::::::::::::::::

First, please **open a terminal** or make sure you are logged into your Vagrant Linux VM.
Assuming that you have a 64-bit system, on Linux, download and install Miniconda 3 with

.. code:: console

    $ wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
    $ bash Miniconda3-latest-Linux-x86_64.sh

On MacOS X, download and install with

.. code:: console

    $ curl https://repo.continuum.io/miniconda/Miniconda3-latest-MacOSX-x86_64.sh -o Miniconda3-latest-MacOSX-x86_64.sh
    $ bash Miniconda3-latest-MacOSX-x86_64.sh

For a 32-bit system, URLs and file names are analogous but without the ``_64``.
When you are asked the question

.. code::

    Do you wish the installer to prepend the Miniconda3 install location to PATH ...? [yes|no]

answer with **yes**.
Along with a minimal Python 3 environment, Miniconda contains the package manager Conda_.
After opening a **new terminal**, you can use the new ``conda`` command to install software packages and create isolated environments to, e.g., use different versions of the same package.
We will later use Conda_ to create an isolated environment with all required software for this tutorial.

Step 2: Preparing a working directory
:::::::::::::::::::::::::::::::::::::

First, **create a new directory** ``snakemake-tutorial`` at a reasonable place and change into that directory in your terminal.
If you use a Vagrant Linux VM from Windows as described above, create that directory under ``/vagrant/``, so that the contents are shared with your host system (you can then edit all files from within Windows with an editor that supports Unix line breaks).
Then, **change to the newly created directory**.
In this directory, we will later create an example workflow that illustrates the Snakemake syntax and execution environment.
First, we download some example data on which the workflow shall be executed:

.. code:: console

    $ wget https://bitbucket.org/snakemake/snakemake-tutorial/get/v3.9.0-1.tar.bz2
    $ tar -xf v3.9.0-1.tar.bz2 --strip 1

This will create a folder ``data`` and a file ``environment.yaml`` in the working directory.

Step 3: Creating an environment with the required software
::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

The ``environment.yaml`` file can be used to install all required software into an isolated Conda environment with the name ``snakemake-tutorial`` via

.. code:: console

    $ conda env create --name snakemake-tutorial --file environment.yaml

Step 4: Activating the environment
::::::::::::::::::::::::::::::::::

To activate the ``snakemake-tutorial`` environment, execute

.. code:: console

    $ source activate snakemake-tutorial

Now you can use the installed tools.
Execute

.. code:: console

    $ snakemake --help

to test this and get information about the command-line interface of Snakemake.
To exit the environment, you can execute

.. code:: console

    $ source deactivate

but **don't do that now**, since we finally want to start working with Snakemake :-).
