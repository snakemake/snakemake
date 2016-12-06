.. _tutorial-setup:

Setup
-----

To go through this tutorial, you need the following software installed:

* Python_ ≥3.3
* Snakemake_ 3.4.2
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

.. code:: bash

    vagrant init hashicorp/precise64
    vagrant up

If you decide to use a 32-bit image, you will need to download the 32-bit version of Miniconda in the next step.
The contents of the ``vagrant-linux`` folder will be shared with the virtual machine that is set up by vagrant.
You can log into the virtual machine via

.. code:: bash

    vagrant ssh

If this command tells you to install an SSH client, you can follow the instructions in this Blogpost_.
Now, you can follow the steps of our tutorial from within your Linux VM.


Step 1: Installing Miniconda 3
::::::::::::::::::::::::::::::

First, please **open a terminal** or make sure you are logged into your Vagrant Linux VM.
Assuming that you have a 64-bit system, on Linux, download and install Miniconda 3 with

.. code:: bash

    wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
    bash Miniconda3-latest-Linux-x86_64.sh

On MacOS X, download and install with

.. code:: bash

    curl https://repo.continuum.io/miniconda/Miniconda3-latest-MacOSX-x86_64.sh -o Miniconda3-latest-MacOSX-x86_64.sh
    bash Miniconda3-latest-MacOSX-x86_64.sh

For a 32-bit system, URLs and file names are analogous but without the ``_64``.
When you are asked the question

.. code::

    Do you wish the installer to prepend the Miniconda3 install location to PATH ...? [yes|no]

answer with **yes**.
Along with a minimal Python 3 environment, Miniconda contains the package manager Conda_.
After opening a **new terminal**, you can use the new ``conda`` command to install software packages and create isolated environments to, e.g., use different versions of the same package.
We will later use Conda_ to create an isolated enviroment with all required software for this tutorial.

Step 2: Preparing a working directory
:::::::::::::::::::::::::::::::::::::

First, **create a new directory** ``snakemake-tutorial`` at a reasonable place and **change into that directory** in your terminal.
If you use a Vagrant Linux VM from Windows as described above, create the directory under ``/vagrant/``, so that the contents are shared with your host system (you can then edit all files from within Windows with an editor that supports Unix line breaks).
In this directory, we will later create an example workflow that illustrates the Snakemake syntax and execution environment.
First, we download some example data on which the workflow shall be executed:

.. code:: bash

    wget https://bitbucket.org/snakemake/snakemake/downloads/snakemake-tutorial-data.tar.gz
    tar -xf snakemake-tutorial-data.tar.gz

This will create a ``data`` folder and a ``requirements.txt`` file in the working directory.

Step 3: Creating an environment with the required software
::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

The ``requirements.txt`` file can be used to install all required software into an isolated conda environment with the name ``snakemake-tutorial`` via

.. code:: bash

    conda create -n snakemake-tutorial -c bioconda --file requirements.txt

Note that the arguments after the ``-c`` flags define software channels that shall be used in addition to the main ``conda`` repository.
Here, we use the Bioconda_ channel, which contains a growing collection of bioinformatics software packaged for Conda.

Step 4: Activating the environment
::::::::::::::::::::::::::::::::::

To activate the ``snakemake-tutorial`` enviroment, execute

.. code:: bash

    source activate snakemake-tutorial

Now you can use the installed tools.
Execute

.. code:: bash

    snakemake --help

to test this and get information about the command-line interface of Snakemake.
To exit the environment, you can execute

.. code:: bash

    source deactivate

but **don't do that now**, since we finally want to start working with Snakemake :-).
