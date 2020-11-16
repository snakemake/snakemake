

================================
Distribution and Reproducibility
================================

It is recommended to store each workflow in a dedicated git repository of the
following structure:

.. code-block:: none

    ├── .gitignore
    ├── README.md
    ├── LICENSE.md
    ├── workflow
    │   ├── rules
    |   │   ├── module1.smk
    |   │   └── module2.smk
    │   ├── envs
    |   │   ├── tool1.yaml
    |   │   └── tool2.yaml
    │   ├── scripts
    |   │   ├── script1.py
    |   │   └── script2.R
    │   ├── notebooks
    |   │   ├── notebook1.py.ipynb
    |   │   └── notebook2.r.ipynb
    │   ├── report
    |   │   ├── plot1.rst
    |   │   └── plot2.rst
    |   └── Snakefile
    ├── config
    │   ├── config.yaml
    │   └── some-sheet.tsv
    ├── results
    └── resources

In other words, the workflow code goes into a subfolder ``workflow``, while the configuration is stored in a subfolder ``config``. 
Inside of the ``workflow`` subfolder, the central ``Snakefile`` marks the entrypoint of the workflow (it will be automatically discovered when running snakemake from the root of above structure. 
In addition to the central ``Snakefile``, rules can be stored in a modular way, using the optional subfolder ``workflow/rules``. Such modules should end with ``.smk`` the recommended file extension of Snakemake.
Further, :ref:`scripts <snakefiles-external_scripts>` should be stored in a subfolder ``workflow/scripts`` and notebooks in a subfolder ``workflow/notebooks``.
Conda environments (see :ref:`integrated_package_management`) should be stored in a subfolder ``workflow/envs`` (make sure to keep them as finegrained as possible to improve transparency and maintainability).
Finally, :ref:`report caption files <snakefiles-reports>` should be stored in ``workflow/report``.
All output files generated in the workflow should be stored under ``results``, unless they are rather retrieved resources, in which case they should be stored under ``resources``. The latter subfolder may also contain small resources that shall be delivered along with the workflow via git (although it might be tempting, please refrain from trying to generate output file paths with string concatenation of a central ``outdir`` variable or so, as this hampers readability).

Then, a workflow can be deployed to a new system via the following steps

.. code-block:: python

    # clone workflow into working directory
    git clone https://github.com/user/myworkflow.git path/to/workdir
    cd path/to/workdir

    # edit config and workflow as needed
    vim config/config.yaml

    # execute workflow, deploy software dependencies via conda
    snakemake -n --use-conda

Importantly, git branching and pull requests can be used to modify and possibly re-integrate workflows.
A `cookiecutter <https://github.com/audreyr/cookiecutter>`_ template for creating this structure can be found `here <https://github.com/snakemake-workflows/cookiecutter-snakemake-workflow>`_.
Given that cookiecutter is installed, you can use it via:

.. code-block:: bash

    cookiecutter gh:snakemake-workflows/cookiecutter-snakemake-workflow

Visit the `Snakemake Workflows Project <https://github.com/snakemake-workflows/docs>`_ for best-practice workflows.

----------------------------------
Uploading workflows to WorkflowHub
----------------------------------

In order to share a workflow with the scientific community it is advised to upload the repository to `WorkflowHub <https://workflowhub.eu/>`_, where each submission will be automatically parsed and encapsulated into a `Research Object Crate <https://w3id.org/ro/crate>`_. That way a *snakemake* workflow is annotated with proper metatada and thus complies with the `FAIR <https://en.wikipedia.org/wiki/FAIR_data>`_ principles of scientific data.

To adhere to the high WorkflowHub standards of scientific workflows the recommended *snakemake* repository structure presented above needs to be extended by the following elements:

- Code of Conduct
- Contribution instructions
- Workflow rule graph
- Workflow documentation
- Test directory

A code of conduct for the repository developers as well as instruction on how to contribute to the project should be placed in the top-level files: ``CODE_OF_CONDUCT.md`` and ``CONTRIBUTING.md``, respectively. Each *snakemake* workflow repository needs to contain an SVG-formatted rule graph placed in a subdirectory ``images/rulegraph.svg``. Additionally, the workflow should be annotated with a technical documentation of all of its subsequent steps, described in ``workflow/documentation.md``. Finally, the repository should contain a ``.tests`` directory with two subdirectories: ``.tests/integration`` and ``.tests/unit``. The former has to contain all the input data, configuration specifications and shell commands required to run an integration test of the whole workflow. The latter shall contain subdirectories dedicated to testing each of the separate workflow steps independently. To simplify the testing procedure *snakemake* can automatically generate unit tests from a successful workflow execution (see :ref:`snakefiles-testing`).

Therefore, the repository structure should comply with:

.. code-block:: none

    ├── .gitignore
    ├── README.md
    ├── LICENSE.md
    ├── CODE_OF_CONDUCT.md
    ├── CONTRIBUTING.md
    ├── .tests
    │   ├── integration
    │   └── unit
    ├── images
    │   └── rulegraph.svg
    ├── workflow
    │   ├── rules
    |   │   ├── module1.smk
    |   │   └── module2.smk
    │   ├── envs
    |   │   ├── tool1.yaml
    |   │   └── tool2.yaml
    │   ├── scripts
    |   │   ├── script1.py
    |   │   └── script2.R
    │   ├── notebooks
    |   │   ├── notebook1.py.ipynb
    |   │   └── notebook2.r.ipynb
    │   ├── report
    |   │   ├── plot1.rst
    |   │   └── plot2.rst
    │   ├── Snakefile
    |   └── documentation.md
    ├── config
    │   ├── config.yaml
    │   └── some-sheet.tsv
    ├── results
    └── resources


.. _integrated_package_management:

-----------------------------
Integrated Package Management
-----------------------------

With Snakemake 3.9.0 it is possible to define isolated software environments per rule.
Upon execution of a workflow, the `Conda package manager <https://conda.pydata.org>`_ is used to obtain and deploy the defined software packages in the specified versions. Packages will be installed into your working directory, without requiring any admin/root priviledges.
Given that conda is available on your system (see `Miniconda <https://conda.pydata.org/miniconda.html>`_), to use the Conda integration, add the ``--use-conda`` flag to your workflow execution command, e.g. ``snakemake --cores 8 --use-conda``.
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

with the following `environment definition <https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#create-env-file-manually>`_:


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

Conda deployment also works well for offline or air-gapped environments. Running ``snakemake --use-conda --conda-create-envs-only`` will only install the required conda environments without running the full workflow. Subsequent runs with ``--use-conda`` will make use of the local environments without requiring internet access.

.. _singularity:

--------------------------
Running jobs in containers
--------------------------

As an alternative to using Conda (see above), it is possible to define, for each rule, a (docker) container to use, e.g.,

.. code-block:: python

    rule NAME:
        input:
            "table.txt"
        output:
            "plots/myplot.pdf"
        container:
            "docker://joseespinosa/docker-r-ggplot2"
        script:
            "scripts/plot-stuff.R"

When executing Snakemake with

.. code-block:: bash

    snakemake --use-singularity

it will execute the job within a container that is spawned from the given image.
Allowed image urls entail everything supported by singularity (e.g., ``shub://`` and ``docker://``).
However, ``docker://`` is preferred, as other container runtimes will be supported in the future (e.g. podman).

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

    container: "docker://continuumio/miniconda3:4.4.10"

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

-------------------------
Using environment modules
-------------------------

In high performace cluster systems (HPC), it can be preferable to use environment modules for deployment of optimized versions of certain standard tools.
Snakemake allows to define environment modules per rule:

.. code-block:: python

    rule bwa:
        input:
            "genome.fa"
            "reads.fq"
        output:
            "mapped.bam"
        conda:
            "envs/bwa.yaml"
        envmodules:
            "bio/bwa/0.7.9",
            "bio/samtools/1.9"
        shell:
            "bwa mem {input} | samtools view -Sbh - > {output}"

Here, when Snakemake is executed with ``snakemake --use-envmodules``, it will load the defined modules in the given order, instead of using the also defined conda environment.
Note that although not mandatory, one should always provide either a conda environment or a container (see above), along with environment module definitions.
The reason is that environment modules are often highly platform specific, and cannot be assumed to be available somewhere else, thereby limiting reproducibility.
By definition an equivalent conda environment or container as a fallback, people outside of the HPC system where the workflow has been designed can still execute it, e.g. by running ``snakemake --use-conda`` instead of ``snakemake --use-envmodules``.

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
