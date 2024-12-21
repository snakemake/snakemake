.. _Mamba: https://github.com/mamba-org/mamba

.. _distribution_and_reproducibility:

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
This main structure and the recommendations below are implemented in `this Snakemake workflow template <https://github.com/snakemake-workflows/snakemake-workflow-template>`_ that you can use to `create your own workflow repository with a single click on "Use this template" <https://github.com/snakemake-workflows/snakemake-workflow-template/generate>`_.
In addition to the central ``Snakefile``, rules can be stored in a modular way, using the optional subfolder ``workflow/rules``.
Such modules should end with ``.smk``, the recommended file extension of Snakemake.
Further, :ref:`scripts <snakefiles-external_scripts>` should be stored in a subfolder ``workflow/scripts`` and notebooks in a subfolder ``workflow/notebooks``.
Conda environments (see :ref:`integrated_package_management`) should be stored in a subfolder ``workflow/envs`` (make sure to keep them as finegrained as possible to improve transparency and maintainability).
Finally, :ref:`report caption files <snakefiles-reports>` should be stored in ``workflow/report``.
All output files generated in the workflow should be stored under ``results``, unless they are rather retrieved resources, in which case they should be stored under ``resources``.
The latter subfolder may also contain small resources that shall be delivered along with the workflow via git (although it might be tempting, please refrain from trying to generate output file paths with string concatenation of a central ``outdir`` variable or so, as this hampers readability).

Workflows set up in above structure can be easily used and combined via :ref:`the Snakemake module system <use_with_modules>`.
Such deployment can even be automated via  `Snakedeploy <https://snakedeploy.readthedocs.io>`_.
Moreover, by publishing a workflow on `Github <https://github.com>`_ and following a set of additional `rules <https://snakemake.github.io/snakemake-workflow-catalog/?rules=true>`_ the workflow will be automatically included in the `Snakemake workflow catalog <https://snakemake.github.io/snakemake-workflow-catalog>`_, thereby easing discovery and even automating its usage documentation.
For an example of such automated documentation, see `here <https://snakemake.github.io/snakemake-workflow-catalog/?usage=snakemake-workflows%2Fdna-seq-varlociraptor>`_.

Visit the `Snakemake Workflows Project <https://github.com/snakemake-workflows/docs>`_ for more best-practice workflows.

.. _use_with_modules:

-----------------------------------------
Using and combining pre-exising workflows
-----------------------------------------

Via the :ref:`module/use <snakefiles-modules>` system introduced with Snakemake 6.0, it is very easy to deploy existing workflows for new projects.
This ranges from the simple application to new data to the complex combination of several complementary workflows in order to perform an integrated analysis over multiple data types.

Consider the following example:

.. code-block:: python

    from snakemake.utils import min_version
    min_version("6.0")

    configfile: "config/config.yaml"

    module dna_seq:
        snakefile:
            # here, it is also possible to provide a plain raw URL like "https://github.com/snakemake-workflows/dna-seq-gatk-variant-calling/raw/v2.0.1/workflow/Snakefile"
            github("snakemake-workflows/dna-seq-gatk-variant-calling", path="workflow/Snakefile", tag="v2.0.1")
        config:
            config

    use rule * from dna_seq

First, we load a local configuration file.
Next, we define the module ``dna_seq`` to be loaded from the URL ``https://github.com/snakemake-workflows/dna-seq-gatk-variant-calling/raw/v2.0.1/workflow/Snakefile``, while using the contents of the local configuration file.
Note that it is possible to either specify the full URL pointing to the raw Snakefile as a string or to use the github marker as done here.
With the latter, Snakemake can however cache the used source files persistently (if a tag is given), such that they don't have to be downloaded on each invocation.
Finally we declare all rules of the dna_seq module to be used.

This kind of deployment is equivalent to just cloning the original repository and modifying the configuration in it.
However, the advantage here is that we are (a) able to easily extend of modify the workflow, while making the changes transparent, and (b) we can store this workflow in a separate (e.g. private) git repository, along with for example configuration and meta data, without the need to duplicate the workflow code.
Finally, we are always able to later combine another module into the current workflow, e.g. when further kinds of analyses are needed.
The ability to modify rules upon using them (see :ref:`snakefiles-modules`) allows for arbitrary rewiring and configuration of the combined modules.

For example, we can easily add another rule to extend the given workflow:

.. code-block:: python

    from snakemake.utils import min_version
    min_version("6.0")

    configfile: "config/config.yaml"

    module dna_seq:
        snakefile:
            # here, it is also possible to provide a plain raw URL like "https://github.com/snakemake-workflows/dna-seq-gatk-variant-calling/raw/v2.0.1/workflow/Snakefile"
            github("snakemake-workflows/dna-seq-gatk-variant-calling", path="workflow/Snakefile", tag="v2.0.1")
        config: config

    use rule * from dna_seq as dna_seq_*

    # easily extend the workflow
    rule plot_vafs:
        input:
            "filtered/all.vcf.gz"
        output:
            "results/plots/vafs.svg"
        notebook:
            "notebooks/plot-vafs.py.ipynb"

    # Define a new default target that collects both the targets from the dna_seq module as well as
    # the new plot.
    rule all:
        input:
            rules.dna_seq_all.input,
            "results/plots/vafs.svg",
        default_target: True

Above, we have added a prefix to all rule names of the dna_seq module, such that there is no name clash with the added rules (``as dna_seq_*`` in the ``use rule`` statement).
In addition, we have added a new rule ``all``, defining the default target in case the workflow is executed (as usually) without any specific target files or rule.
The new target rule collects both all input files of the rule ``all`` from the dna_seq workflow, as well as additionally collecting the new plot.

It is possible to further extend the workflow with other modules, thereby generating an integrative analysis.
Here, let us assume that we want to conduct another kind of analysis, say RNA-seq, using a different external workflow.
We can extend above example in the following way:

.. code-block:: python

    from snakemake.utils import min_version
    min_version("6.0")

    configfile: "config/config.yaml"

    module dna_seq:
        snakefile:
            github("snakemake-workflows/dna-seq-gatk-variant-calling", path="workflow/Snakefile", tag="v2.0.1")
        config: config["dna-seq"]
        prefix: "dna-seq"

    use rule * from dna_seq as dna_seq_*

    rule plot_vafs:
        input:
            "filtered/all.vcf.gz"
        output:
            "results/plots/vafs.svg"
        notebook:
            "notebooks/plot-vafs.py.ipynb"

    module rna_seq:
        snakefile:
            github("snakemake-workflows/rna-seq-kallisto-sleuth", path="workflow/Snakefile", tag="v2.0.1")
        config: config["rna-seq"]
        prefix: "rna-seq"

    use rule * from rna_seq as rna_seq_*


    # Define a new default target that collects all the targets from the dna_seq and rna_seq module.
    rule all:
        input:
            rules.dna_seq_all.input,
            rules.rna_seq_all.input,
        default_target: True

Above, several things have changed.

* First, we have added another module ``rna_seq``.
* Second, we have added a prefix to all non-absolute input and output file names of both modules (``prefix: "dna-seq"`` and ``prefix: "rna-seq"``) in order to avoid file name clashes.
* Third, we have added a default target rule that collects both the default targets from the module ``dna_seq`` as well as the module ``rna_seq``.
* Finally, we provide the config of the two modules via two separate sections in the common config file (``config["dna-seq"]`` and ``config["rna-seq"]``).

----------------------------------
Uploading workflows to WorkflowHub
----------------------------------

In order to share a workflow with the scientific community it is advised to upload the repository to `WorkflowHub <https://workflowhub.eu/>`_, where each submission will be automatically parsed and encapsulated into a `Research Object Crate <https://w3id.org/ro/crate>`_. That way a *snakemake* workflow is annotated with proper metadata and thus complies with the `FAIR <https://en.wikipedia.org/wiki/FAIR_data>`_ principles of scientific data.

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

It is possible (and highly encouraged, see :ref:`snakefiles-best_practices`) to define isolated software environments per rule.
Upon execution of a workflow, the `Conda package manager <https://conda.pydata.org>`_ is used to obtain and deploy the defined software packages in the specified versions. Packages will be installed into your working directory, without requiring any admin/root privileges.
Given that conda is available on your system (see `Miniconda <https://conda.pydata.org/miniconda.html>`_), to use the Conda integration, add the ``--software-deployment-method conda`` option (``--sdm`` for short) to your workflow execution command, e.g. ``snakemake --cores 8 --sdm conda``.
When ``--software-deployment-method conda`` (``--sdm`` for short) is activated, Snakemake will automatically create software environments for any used wrapper (see :ref:`snakefiles-wrappers`).
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

Please note that in the environment definition, conda determines the priority of channels depending on their order of appearance in the channels list. For instance, the channel that comes first in the list gets the highest priority. Default packages defined in the user configuration of conda (`.condarc`)) are ignored by Snakemake.

The path to the environment definition is interpreted as **relative to the Snakefile that contains the rule** (unless it is an absolute path, which is discouraged).

Instead of using a concrete path, it is also possible to provide a path containing wildcards (which must also occur in the output files of the rule), analogous to the specification of input files.

In addition, it is possible to use a callable which returns a ``str`` value.
The signature of the callable has to be ``callable(wildcards [, params] [, input])`` (``params`` and ``input`` are optional parameters).

Note that the use of distinct conda environments for different jobs from the same rule is currently not properly displayed in the generated reports.
At the moment, only a single, random conda environment is shown.

.. note::

   Note that conda environments are only used with ``shell``, ``script``, ``notebook`` and the ``wrapper`` directive, not the ``run`` directive.
   The reason is that the ``run`` directive has access to the rest of the Snakefile (e.g. globally defined variables) and therefore must be executed in the same process as Snakemake itself. If used with ``notebook`` directive, the associated conda environment should have package ``jupyter`` installed (this package contains dependencies required to execute the notebook).

   Further, note that search path modifying environment variables like ``R_LIBS`` and ``PYTHONPATH`` can interfere with your conda environments.
   Therefore, Snakemake automatically deactivates them for a job when a conda environment definition is used.
   If you know what you are doing, in order to deactivate this behavior, you can use the flag ``--conda-not-block-search-path-envvars``.

Snakemake will store the environment persistently in ``.snakemake/conda/$hash`` with ``$hash`` being the MD5 hash of the environment definition file content. This way, updates to the environment definition are automatically detected.
Note that you need to clean up environments manually for now. However, in many cases they are lightweight and consist of symlinks to your central conda installation.

Conda deployment also works well for offline or air-gapped environments. Running ``snakemake --sdm conda --conda-create-envs-only`` will only install the required conda environments without running the full workflow. Subsequent runs with ``--sdm conda`` will make use of the local environments without requiring internet access.

Freezing environments to exactly pinned packages
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If Snakemake finds a special file ending on ``<platform>.pin.txt`` next to a conda environment file (with ``<platform>`` being the current platform, e.g. ``linux-64``), it will try to use the contents of that file to determine the conda packages to deploy.
The file is expected to contain conda's `explicit specification file format <https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#building-identical-conda-environments>`_.
Snakemake will first try to deploy the environment using that file, and only if that fails it will use the regular environment file.

This enables to freeze an environment to a certain state, and will ensure that people using a workflow will get exactly the same environments down to the individual package builds, which is in fact very similar to providing the environment encapsulated in a container image.
Generating such pin files for conda environments can be automatically done using `Snakedeploy <https://snakedeploy.readthedocs.io>`_.
Let ``envs/ggplot.yaml`` be the conda environment file used in the example above.
Then, the pinning can be generated with

.. code-block:: bash

    snakedeploy pin-conda-envs envs/ggplot.yaml

Multiple paths to environments can be provided at the same time; also see ``snakedeploy pin-conda-envs --help``.

Of course, it is **important to update the pinnings** whenever the original environment is modified, such that they do not diverge.

Updating environments
~~~~~~~~~~~~~~~~~~~~~

When a workflow contains many conda environments, it can be helpful to automatically update them to the latest versions of all packages.
This can be done automatically via `Snakedeploy <https://snakedeploy.readthedocs.io>`_:

.. code-block:: bash

    snakedeploy update-conda-envs envs/ggplot.yaml

Multiple paths to environments can be provided at the same time; also see ``snakedeploy update-conda-envs --help``.


Providing post-deployment scripts
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

From Snakemake 6.14 onwards post-deployment shell-scripts can be provided to perform additional adjustments of a conda environment.
This might be helpful in case a conda package is missing components or requires further configuration for execution.
Post-deployment scripts must be placed next to their corresponding environment-file and require the suffix ``.post-deploy.sh``, e.g.:

.. code-block:: python

    rule NAME:
        input:
            "seqs.fastq"
        output:
            "results.tsv"
        conda:
            "envs/interproscan.yaml"
        shell:
            "interproscan.sh -i {input} -f tsv -o {output}"

.. code-block:: none

    ├── Snakefile
    └── envs
        ├── interproscan.yaml
        └── interproscan.post-deploy.sh

The path of the conda environment can be accessed within the script via ``$CONDA_PREFIX``.
Importantly, if the script relies on certain shell specific syntax, (e.g. `set -o pipefail` for bash), make sure to add a matching shebang to the script, e.g.:

.. code-block:: bash

    #!env bash
    set -o pipefail
    # ...

If no shebang line like above (``#!env bash``) is provided, the script will be executed with the ``sh`` command.

.. _conda_named_env:

-----------------------------------------
Using already existing conda environments
-----------------------------------------

Sometimes it can be handy to refer to an already existing conda environment from a rule, instead of defining a new one from scratch.
Importantly, one should be aware that this can **hamper reproducibility**, because the workflow then relies on this environment to be present
**in exactly the same way** on any new system where the workflow is executed. Essentially, you will have to take care of this manually in such a case.
Therefore, the approach using environment definition files described above is highly recommended and preferred.
Referring to an existing environment can however be useful during development, e.g. when a certain software package is developed in parallel to a workflow that uses it.

It is possible to refer to a named environment:

.. code-block:: python

    rule NAME:
        input:
            "table.txt"
        output:
            "plots/myplot.pdf"
        conda:
            "some-env-name"
        script:
            "scripts/plot-stuff.R"

Or alternatively to the filesystem path of an environment:

.. code-block:: python

    rule NAME:
        input:
            "table.txt"
        output:
            "plots/myplot.pdf"
        conda:
            "/home/johannes/miniforge/envs/some-env-name"
        script:
            "scripts/plot-stuff.R"

For any such rules, Snakemake will just activate the given environment, instead of automatically deploying anything.
Instead of using a concrete name or path, it is also possible to provide one containing wildcards (which must also occur in the output files of the rule), analogous to the specification of input files.
Finally it is also possible to use a callable which returns a ``str`` value and takes ``wildcards`` as single argument, similar to input functions.

Note that Snakemake distinguishes file based environments from named ones as follows:
if the given specification ends on ``.yaml`` or ``.yml``, Snakemake assumes it to be a path to an environment definition file;
otherwise, it assumes the given specification to point to an existing environment.


.. _apptainer:

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
            "docker://joseespinosa/docker-r-ggplot2:1.0"
        script:
            "scripts/plot-stuff.R"

When executing Snakemake with

.. code-block:: bash

    snakemake --software-deployment-method apptainer
    # or the shorthand version
    snakemake --sdm apptainer

it will execute the job within a container that is spawned from the given image.
Allowed image urls entail everything supported by apptainer (e.g., ``shub://`` and ``docker://``).
However, ``docker://`` is preferred, as other container runtimes will be supported in the future (e.g. podman).

Defining global container images
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. note::

   Note that apptainer integration is only used with ``shell``, ``script`` and the ``wrapper`` directive, not the ``run`` directive.
   The reason is that the ``run`` directive has access to the rest of the Snakefile (e.g. globally defined variables) and therefore must be executed in the same process as Snakemake itself.

A global definition of a container image can be given:

.. code-block:: python

    container: "docker://joseespinosa/docker-r-ggplot2:1.0"

    rule NAME:
        ...

In this case all jobs will be executed in a container. You can disable execution in container
by setting the container directive of the rule to ``None``.

.. code-block:: python

    container: "docker://joseespinosa/docker-r-ggplot2:1.0"

    rule NAME:
        container: None


Handling shell executable
~~~~~~~~~~~~~~~~~~~~~~~~~

Snakemake executes rules using the bash shell by default.
If your container image does not contain bash, you can specify a different shell executable, see :ref:`shell_settings`.

-----------------------------------------
Containerization of Conda based workflows
-----------------------------------------
While :ref:`integrated_package_management` provides control over the used software in exactly
the desired versions, it does not control the underlying operating system.
However, given a workflow with conda environments for each rule, Snakemake can automatically
generate a container image specification (in the form of a ``Dockerfile``) that contains
all required environments via the flag --containerize:

.. code-block:: bash

    snakemake --containerize > Dockerfile

The container image specification generated by Snakemake aims to be transparent and readable, e.g. by displaying each contained environment in a human readable way.
Via the special directive ``containerized`` this container image can be used in the workflow (both globally or per rule) such that no further conda package downloads are necessary, for example:

.. code-block:: python

    containerized: "docker://username/myworkflow:1.0.0"

    rule NAME:
        input:
            "table.txt"
        output:
            "plots/myplot.pdf"
        conda:
            "envs/ggplot.yaml"
        script:
            "scripts/plot-stuff.R"

Using the containerization of Snakemake has three advantages over manually crafting a container image for a workflow:

1. A workflow with conda environment definitions is much more transparent to the reader than a black box container image, as each rule directly shows which software stack is used. Containerization just persistently projects those environments into a container image.
2. It remains possible to run the workflow without containers, just via the conda environments.
3. During development, testing can first happen without the container and just on the conda environments. When releasing a production version of the workflow the image can be uploaded just once and for future stable releases, thereby limiting the overhead created in container registries.

--------------------------------------------------------------
Ad-hoc combination of Conda package management with containers
--------------------------------------------------------------

While :ref:`integrated_package_management` provides control over the used software in exactly
the desired versions, it does not control the underlying operating system.
Here, it becomes handy that Snakemake >=4.8.0 allows to combine Conda-based package management
with :ref:`apptainer`.
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

    snakemake --software-deployment-method conda apptainer
    # or the shorthand version
    snakemake --sdm conda apptainer

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

In high performance cluster systems (HPC), it can be preferable to use environment modules for deployment of optimized versions of certain standard tools.
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
By definition an equivalent conda environment or container as a fallback, people outside of the HPC system where the workflow has been designed can still execute it, e.g. by running ``snakemake --software-deployment-method conda`` instead of ``snakemake --use-envmodules``.

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

.. _global-workflow-dependencies:

----------------------------
Global workflow dependencies
----------------------------

Often, your workflow will depend on some additional packages that need to be present 
along with Snakemake in order to handle actions before any rule is executed.
Classical examples for this are `pandas <https://pandas.pydata.org/>`_, 
`pep <https://pep.databio.org>`_ (also see :ref:`snakefiles-peps`) and 
:ref:`storage plugins <storage-support>`.

Snakemake allows to define such global dependencies using a global ``conda`` directive
that should occur at the beginning of your workflow, before you import or use any of
those additional packages::

    conda:
        "envs/global.yaml"

With ``envs/global.yaml`` containing e.g.::

    channels:
      - conda-forge
      - bioconda
      - nodefaults
    dependencies:
      - pandas=1.0.3
      - snakemake-storage-plugin-s3
    
Under the hood, this is implemented using `conda-inject <https://github.com/koesterlab/conda-inject>`_, 
which modifies the python searchpath and the PATH variable on the fly during execution,
pointing to additional environments that do not alter the environment in which Snakemake
has been installed.

This mechanism requires that you use Mamba_ or Conda and activate conda-based software deployment via::

    --software-deployment-method conda
    # or the shorthand version
    --sdm conda
