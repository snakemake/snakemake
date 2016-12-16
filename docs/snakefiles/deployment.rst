

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
When ``--use-conda`` is activated, Snakemake will automatically create software environments for any used wrapper (see above).
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
Note that you need to clean up environments manually for now. However, they are lightweight and consist only of symlinks to your central conda installation.
