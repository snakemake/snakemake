
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
.. _Miniforge: https://github.com/conda-forge/miniforge
.. _Mamba: https://github.com/mamba-org/mamba
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

However, don't install any of these this manually now, we guide you through better ways below.

.. _tutorial-free-on-gitpod:

Run tutorial for free in the cloud via Gitpod
:::::::::::::::::::::::::::::::::::::::::::::

.. note::

    A common thing to happen while using the development environment in GitPod is to hit ``Ctrl-s`` while in the terminal window, because you wanted to save a file in the editor window.
    This will freeze up you terminal.
    To get it back, make sure you selected the terminal window by clicking on it and then hit ``Ctrl-q``.

The easiest way to run this tutorial is to use Gitpod, which enables performing the exercises via your browser---including all required software, for free and in the cloud.
In order to do this, simply open the predefined `snakemake-tutorial GitPod workspace <https://gitpod.io/#https://github.com/snakemake/snakemake-tutorial-data>`_ in your browser.
GitPod provides you with a `Theia development environment <https://theia-ide.org/docs>`_, which you can learn about in the linked documentation.
Once you have a basic understanding of this environment, you can go on directly with :ref:`tutorial-basics`.

Running the tutorial on your local machine
::::::::::::::::::::::::::::::::::::::::::

If you prefer to run the tutorial on your local machine, please follow the steps below.

The easiest way to set these prerequisites up, is to use the Miniforge_ Python 3 distribution
(Miniforge_ is a Conda based distribution like Miniconda_, which however uses Mamba_ a fast and more robust replacement for the Conda_ package manager).
The tutorial assumes that you are using either Linux or MacOS X.
Both Snakemake and Miniforge_ work also under Windows, but the Windows shell is too different to be able to provide generic examples.

**Currently, the setup currently only works for Intel based machines (x86_64), not ARM based machines like the new Apple M1/2/3 architecture.**
This will change in the coming months. In the meantime, if you are on an ARM based Mac, you can use Rosetta to emulate an intel architecture.
Otherwise, you can simply use the Gitpod approach outlined above.

Setup on Windows
::::::::::::::::

If you already use Linux or MacOS X, go on with **Step 1**.

Windows Subsystem for Linux
"""""""""""""""""""""""""""

If you use Windows 10, you can set up the Windows Subsystem for Linux (`WSL`_) to natively run linux applications.
Install the WSL following the instructions in the `WSL Documentation`_. You can chose any Linux distribution available for the WSL, but the most popular and accessible one is Ubuntu.
Start the WSL and set up your account; now, you can follow the steps of our tutorial from within your Linux environment in the WSL.

Vagrant virtual machine
"""""""""""""""""""""""

If you are using a version of Windows older than 10 or if you do not wish to install the WSL, you can instead setup a Linux virtual machine (VM) with Vagrant_.
First, install Vagrant following the installation instructions in the `Vagrant Documentation`_.
Then, create a new directory you want to share with your Linux VM, for example, create a folder named ``vagrant-linux`` somewhere.
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


Step 1: Installing Miniforge
:::::::::::::::::::::::::::::

First, please **open a terminal** or make sure you are logged into your Vagrant Linux VM.
Assuming that you have a 64-bit system, on Linux, download and install Miniconda 3 with

.. code:: console

    $ curl -L https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh -o Miniforge3-Linux-x86_64.sh
    $ bash Miniforge3-Linux-x86_64.sh

On MacOS with x86_64 architecture, download and install with

.. code:: console

    $ curl -L https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-MacOSX-x86_64.sh -o Miniforge3-MacOSX-x86_64.sh
    $ bash Miniforge3-MacOSX-x86_64.sh

On MacOS with ARM/M1 architecture, download and install with

.. code:: console

    $ curl -L https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-MacOSX-arm64.sh -o Miniforge3-MacOSX-arm64.sh
    $ bash Miniforge3-MacOSX-arm64.sh

When you are asked the question

.. code::

    Do you wish the installer to prepend the install location to PATH ...? [yes|no]

answer with **yes**.
Along with a minimal Python 3 environment, Miniforge contains the package manager Mamba_.
After closing your current terminal and opening a **new terminal**, you can use the new ``conda`` command to install software packages and create isolated environments to, for example, use different versions of the same package.
We will later use Conda_ to create an isolated environment with all the required software for this tutorial.

Step 2: Preparing a working directory
:::::::::::::::::::::::::::::::::::::

First, **create a new directory** ``snakemake-tutorial`` at a **place you can easily remember** and change into that directory in your terminal:

.. code:: console

    $ mkdir snakemake-tutorial
    $ cd snakemake-tutorial

If you use a Vagrant Linux VM from Windows as described above, create that directory under ``/vagrant/``, so that the contents are shared with your host system (you can then edit all files from within Windows with an editor that supports Unix line breaks).
Then, **change to the newly created directory**.
In this directory, we will later create an example workflow that illustrates the Snakemake syntax and execution environment.
First, we download some example data on which the workflow shall be executed:

.. code:: console

    $ curl -L https://api.github.com/repos/snakemake/snakemake-tutorial-data/tarball -o snakemake-tutorial-data.tar.gz

Next we extract the data. On Linux, run

.. code:: console

    $ tar --wildcards -xf snakemake-tutorial-data.tar.gz --strip 1 "*/data" "*/environment.yaml"

On MacOS, run

.. code:: console

    $ tar -xf snakemake-tutorial-data.tar.gz --strip 1 "*/data" "*/environment.yaml"

This will create a folder ``data`` and a file ``environment.yaml`` in the working directory.

Step 3: Creating an environment with the required software
::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

All interactions with Conda package management below can be conducted with either ``conda``, ``mamba`` or ``micromamba``.
For the steps in the :ref:`"advanced" part of the tutorial <tutorial-advanced>`, you have to install ``mamba`` though in case you don't have it.

First, make sure to activate the base environment with

.. code:: console

    $ conda activate base

The ``environment.yaml`` file that you have obtained with the previous step (Step 2) can be used to install all required software into an isolated Conda environment with the name ``snakemake-tutorial`` via

.. code:: console

    $ mamba env create --name snakemake-tutorial --file environment.yaml

If you don't have the Mamba_ command because you used a different conda distribution than Miniforge_, you can also first install Mamba_
(which is a faster and more robust replacement for Conda_) in your base environment with

.. code:: console

    $ conda install -n base -c conda-forge mamba

and then run the `mamba env create` command shown above.

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
