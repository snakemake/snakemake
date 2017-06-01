

====================================
Workflow Distribution and Deployment
====================================

It is recommended to store each workflow in a dedicated git repository of the
following structure:

.. code-block:: none

    ├── .gitignore
    ├── README.md
    ├── LICENSE.md
    ├── config.yaml
    ├── environment.yaml
    ├── scripts
    │   ├── __init__.py
    │   ├── script1.py
    │   └── script2.R
    └── Snakefile

Then, a workflow can be deployed to a new system via the following steps

.. code-block:: python

    # clone workflow into working directory
    git clone https://bitbucket.org/user/myworkflow.git path/to/workdir
    cd path/to/workdir

    # edit config and workflow as needed
    vim config.yaml

    # install dependencies into isolated environment
    conda env create -n myworkflow --file environment.yaml

    # activate environment
    source activate myworkflow

    # execute workflow
    snakemake -n

Importantly, git branching and pull requests can be used to modify and possibly re-integrate workflows.

.. _integrated_package_management:

-----------------------------
Integrated Package Management
-----------------------------

With Snakemake 3.9.0 it is possible to define isolated software environments per rule.
Upon execution of a workflow, the `Conda package manager <http://conda.pydata.org>`_ is used to obtain and deploy the defined software packages in the specified versions. Packages will be installed into your working directory, without requiring any admin/root priviledges.
Given that conda is available on your system (see `Miniconda <http://conda.pydata.org/miniconda.html>`_), to use the Conda integration, add the ``--use-conda`` flag to your workflow execution command, e.g. ``snakemake --cores 8 --use-conda``.
When ``--use-conda`` is activated, Snakemake will automatically create software environments for any used wrapper (see :ref:`snakefiles-wrappers`).
Further, you can manually define environments via the ``conda`` directive, e.g.:

.. code-block:: python

    rule NAME:
        input:
            "table.txt"
        output:
            "plots/myplot.pdf"
        conda:
            "envs/ggplot.yaml"
        script:
            "scripts/plot-stuff.R"

with the following `environment definition <http://conda.pydata.org/docs/using/envs.html#create-environment-file-by-hand>`_:


.. code-block:: yaml

    channels:
     - r
    dependencies:
     - r=3.3.1
     - r-ggplot2=2.1.0

Snakemake will store the environment persistently in ``.snakemake/conda/$hash`` with ``$hash`` being the MD5 hash of the environment definition file content. This way, updates to the environment definition are automatically detected.
Note that you need to clean up environments manually for now. However, in many cases they are lightweight and consist of symlinks to your central conda installation.


--------------------------------------
Sustainable and reproducible archiving
--------------------------------------

With Snakemake 3.10.0 it is possible to archive a workflow into a
`tarball <https://en.wikipedia.org/wiki/Tar_(computing)>`_
(`.tar`, `.tar.gz`, `.tar.bz2`, `.tar.xz`), via

.. code-block:: bash

    snakemake --archive my-workflow.tar.gz

If above layout is followed, this will archive any code and config files that
is under git version control. Further, all input files will be included into the
archive. Finally, the software packages of each defined conda environment are included.
This results in a self-contained workflow archive that can be re-executed on a
vanilla machine that only has Conda and Snakemake installed via

.. code-block:: bash

    tar -xf my-workflow.tar.gz
    snakemake -n

Note that the archive is platform specific. For example, if created on Linux, it will
run on any Linux newer than the minimum version that has been supported by the used
Conda packages at the time of archiving (e.g. CentOS 6).

A useful pattern when publishing data analyses is to create such an archive,
upload it to `Zenodo <https://zenodo.org/>`_ and thereby obtain a
`DOI <https://en.wikipedia.org/wiki/Digital_object_identifier>`_.
Then, the DOI can be cited in manuscripts, and readers are able to download
and reproduce the data analysis at any time in the future.
