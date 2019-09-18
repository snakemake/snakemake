

================================
Distribution and Reproducibility
================================

It is recommended to store each workflow in a dedicated git repository of the
following structure:

.. code-block:: none

    ├── .gitignore
    ├── README.md
    ├── LICENSE.md
    ├── config.yaml
    ├── scripts
    │   ├── script1.py
    │   └── script2.R
    ├── envs
    │   └── myenv.yaml
    └── Snakefile

Then, a workflow can be deployed to a new system via the following steps

.. code-block:: python

    # clone workflow into working directory
    git clone https://bitbucket.org/user/myworkflow.git path/to/workdir
    cd path/to/workdir

    # edit config and workflow as needed
    vim config.yaml

    # execute workflow, deploy software dependencies via conda
    snakemake -n --use-conda

Importantly, git branching and pull requests can be used to modify and possibly re-integrate workflows.
A `cookiecutter <https://github.com/audreyr/cookiecutter>`_ template for creating this structure can be found `here <https://github.com/snakemake-workflows/cookiecutter-snakemake-workflow>`_.
Given that cookiecutter is installed, you can use it via:

.. code-block:: bash

    cookiecutter gh:snakemake-workflows/cookiecutter-snakemake-workflow

Visit the `Snakemake Workflows Project <https://github.com/snakemake-workflows/docs>`_ for best-practice workflows.

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

The path to the environment definition is interpreted as **relative to the Snakefile that contains the rule** (unless it is an absolute path, which is discouraged).

.. sidebar:: Note

   Note that conda environments are only used with ``shell``, ``script`` and the ``wrapper`` directive, not the ``run`` directive.
   The reason is that the ``run`` directive has access to the rest of the Snakefile (e.g. globally defined variables) and therefore must be executed in the same process as Snakemake itself.

Snakemake will store the environment persistently in ``.snakemake/conda/$hash`` with ``$hash`` being the MD5 hash of the environment definition file content. This way, updates to the environment definition are automatically detected.
Note that you need to clean up environments manually for now. However, in many cases they are lightweight and consist of symlinks to your central conda installation.


.. _singularity:

--------------------------
Running jobs in containers
--------------------------

As an alternative to using Conda (see above), it is possible to define, for each rule, a docker or singularity container to use, e.g.,

.. code-block:: python

    rule NAME:
        input:
            "table.txt"
        output:
            "plots/myplot.pdf"
        singularity:
            "docker://joseespinosa/docker-r-ggplot2"
        script:
            "scripts/plot-stuff.R"

When executing Snakemake with

.. code-block:: bash

    snakemake --use-singularity

it will execute the job within a singularity container that is spawned from the given image.
Allowed image urls entail everything supported by singularity (e.g., ``shub://`` and ``docker://``).

.. sidebar:: Note

   Note that singularity integration is only used with ``shell``, ``script`` and the ``wrapper`` directive, not the ``run`` directive.
   The reason is that the ``run`` directive has access to the rest of the Snakefile (e.g. globally defined variables) and therefore must be executed in the same process as Snakemake itself.


When ``--use-singularity`` is combined with ``--kubernetes`` (see :ref:`kubernetes`), cloud jobs will be automatically configured to run in priviledged mode, because this is a current requirement of the singularity executable.
Importantly, those privileges won't be shared by the actual code that is executed in the singularity container though.

--------------------------------------------------
Combining Conda package management with containers
--------------------------------------------------

While :ref:`integrated_package_management` provides control over the used software in exactly
the desired versions, it does not control the underlying operating system.
Here, it becomes handy that Snakemake >=4.8.0 allows to combine Conda-based package management
with :ref:`singularity`.
For example, you can write

.. code-block:: python

    singularity: "docker://continuumio/miniconda3:4.4.10"

    rule NAME:
        input:
            "table.txt"
        output:
            "plots/myplot.pdf"
        conda:
            "envs/ggplot.yaml"
        script:
            "scripts/plot-stuff.R"

in other words, a global definition of a container image can be combined with a
per-rule conda directive.
Then, upon invocation with

.. code-block:: bash

    snakemake --use-conda --use-singularity

Snakemake will first pull the defined container image, and then create the requested conda environment from within the container.
The conda environments will still be stored in your working environment, such that they don't have to be recreated unless they have changed.
The hash under which the environments are stored includes the used container image url, such that changes to the container image also lead to new environments to be created.
When a job is executed, Snakemake will first enter the container and then activate the conda environment.

By this, both packages and OS can be easily controlled without the overhead of creating and distributing specialized container images.
Of course, it is also possible (though less common) to define a container image per rule in this scenario.

The user can, upon execution, freely choose the desired level of reproducibility:

* no package management (use whatever is on the system)
* Conda based package management (use versions defined by the workflow developer)
* Conda based package management in containerized OS (use versions and OS defined by the workflow developer)

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
