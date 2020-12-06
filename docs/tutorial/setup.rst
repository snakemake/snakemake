
.. _tutorial-setup:

Setup
-----

.. _Snakemake: https://snakemake.readthedocs.io
.. _Snakemake homepage: https://snakemake.readthedocs.io
.. _GNU Make: https://www.gnu.org/software/make
.. _Python: https://www.python.org
.. _BWA: http://bio-bwa.sourceforge.net
.. _SAMtools: https://www.htslib.org
.. _BCFtools: https://www.htslib.org
.. _Pandas: https://pandas.pydata.org
.. _Miniconda: https://conda.pydata.org/miniconda.html
.. _Conda: https://conda.pydata.org
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

Requirements
::::::::::::

To go through this tutorial, you need the following software installed:

* Python_ ≥3.5
* Snakemake_ ≥5.24.1
* BWA_ 0.7
* SAMtools_ 1.9
* Pysam_ 0.15
* BCFtools_ 1.9
* Graphviz_ 2.42
* Jinja2_ 2.11
* NetworkX_ 2.5
* Matplotlib_ 3.3

The easiest way to setup these prerequisites is to use the Miniconda_ Python 3 distribution.
The tutorial assumes that you are using either Linux or MacOS X.
Both Snakemake and Miniconda work also under Windows, but the Windows shell is too different to be able to provide generic examples.

Setup on Windows
::::::::::::::::

If you already use Linux or MacOS X, go on with **Step 1**.

Windows Subsystem for Linux
"""""""""""""""""""""""""""

If you use Windows 10, you can setup the Windows Subsystem for Linux (`WSL`_) to natively run linux applications.
Install the WSL following the instructions in the `WSL Documentation`_. You can chose any Linux distribution available for the WSL, but the most popular and accessible one is Ubuntu.
Start the WSL and setup your account; now, you can follow the steps of our tutorial from within your Linux environment in the WSL.

Vagrant virtual machine
"""""""""""""""""""""""

If you are using a version of Windows older than 10 or if you do not wish to install the WSL, you can instead setup a Linux virtual machine (VM) with Vagrant_.
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

    $ wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
    $ bash Miniconda3-latest-Linux-x86_64.sh

On MacOS X, download and install with

.. code:: console

    $ curl https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh -o Miniconda3-latest-MacOSX-x86_64.sh
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

First, **create a new directory** ``snakemake-tutorial`` at a **reasonable place** and change into that directory in your terminal:

.. code:: console

    $ mkdir snakemake-tutorial
    $ cd snakemake-tutorial

If you use a Vagrant Linux VM from Windows as described above, create that directory under ``/vagrant/``, so that the contents are shared with your host system (you can then edit all files from within Windows with an editor that supports Unix line breaks).
Then, **change to the newly created directory**.
In this directory, we will later create an example workflow that illustrates the Snakemake syntax and execution environment.
First, we download some example data on which the workflow shall be executed:

.. code:: console

    $ wget https://github.com/snakemake/snakemake-tutorial-data/archive/v5.24.1.tar.gz
    $ tar --wildcards -xf v5.24.1.tar.gz --strip 1 "*/data" "*/environment.yaml"

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

    $ conda activate snakemake-tutorial

Now you can use the installed tools.
Execute

.. code:: console

    $ snakemake --help

to test this and get information about the command-line interface of Snakemake.
To exit the environment, you can execute

.. code:: console

    $ conda deactivate

but **don't do that now**, since we finally want to start working with Snakemake :-).
