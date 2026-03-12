.. _Snakemake Wrapper Repository: https://snakemake-wrappers.readthedocs.io

.. _snakefiles-modularization:

==============
Modularization
==============

Modularization in Snakemake comes at four different levels.

1. The most fine-grained level are wrappers. They are available and can be published at the `Snakemake Wrapper Repository`_. These wrappers can then be composed and customized according to your needs, by copying skeleton rules into your workflow. In combination with conda integration, wrappers also automatically deploy the needed software dependencies into isolated environments.
2. For larger, reusable parts that shall be integrated into a common workflow, it is recommended to write small Snakefiles and include them into a main Snakefile via the include statement. In such a setup, all rules share a common config file.
3. The third level is provided via the :ref:`module statement <snakefiles-modules>`, which enables arbitrary combination and reuse of rules.


.. _snakefiles-wrappers:

--------
Wrappers
--------

The wrapper directive allows to have reusable wrapper scripts around e.g. command line tools.
In contrast to modularization strategies like ``include`` or subworkflows, the wrapper directive allows to re-wire the DAG of jobs.
For example

.. code-block:: python

    rule samtools_sort:
        input:
            "mapped/{sample}.bam"
        output:
            "mapped/{sample}.sorted.bam"
        params:
            "-m 4G"
        threads: 8
        wrapper:
            "0.0.8/bio/samtools/sort"

.. note::

    It is possible to refer to wildcards and params in the wrapper identifier, e.g. by specifying ``"0.0.8/bio/{params.wrapper}"`` or ``"0.0.8/bio/{wildcards.wrapper}"``.

Refers to the wrapper ``"0.0.8/bio/samtools/sort"`` to create the output from the input.
Snakemake will automatically download the wrapper from the `Snakemake Wrapper Repository`_.
Thereby, ``0.0.8`` can be replaced with the git `version tag <https://github.com/snakemake/snakemake-wrappers/releases>`_ you want to use, or a `commit id <https://github.com/snakemake/snakemake-wrappers/commits>`_.
This ensures reproducibility since changes in the wrapper implementation will only be propagated to your workflow once you update the version tag.
Examples for each wrapper can be found in the READMEs located in the wrapper subdirectories at the `Snakemake Wrapper Repository`_.

Alternatively, for example during development, the wrapper directive can also point to full URLs, including URLs to local files with absolute paths ``file://`` or relative paths ``file:``.
Such a URL will have to point to the folder containing the ``wrapper.*`` and ``environment.yaml`` files.
In the above example, the full GitHub URL could for example be provided with ``wrapper: https://github.com/snakemake/snakemake-wrappers/raw/0.0.8/bio/samtools/sort``.
Note that it needs to point to the ``/raw/`` version of the folder, not the rendered HTML version.

In addition, the `Snakemake Wrapper Repository`_ offers so-called meta-wrappers, which can be used as modules, see :ref:`snakefiles-meta-wrappers`.

The `Snakemake Wrapper Repository`_ is meant as a collaborative project and pull requests are very welcome.


.. _cwl:

--------------------------------------
Common-Workflow-Language (CWL) support
--------------------------------------

With Snakemake 4.8.0, it is possible to refer to `CWL <https://www.commonwl.org/>`__ tool definitions in rules instead of specifying a wrapper or a plain shell command.
A CWL tool definition can be used as follows.

.. code-block:: python

    rule samtools_sort:
        input:
            input="mapped/{sample}.bam"
        output:
            output_name="mapped/{sample}.sorted.bam"
        params:
            threads=lambda wildcards, threads: threads,
            memory="4G"
        threads: 8
        cwl:
            "https://github.com/common-workflow-language/workflows/blob/"
            "fb406c95/tools/samtools-sort.cwl"

.. note::

    It is possible to refer to wildcards and params in the tool definition URL, e.g. by specifying something like ``"https://.../tools/{params.tool}.cwl"`` or ``"https://.../tools/{wildcards.tool}.cwl"``.

It is advisable to use a github URL that includes the commit as above instead of a branch name, in order to ensure reproducible results.
Snakemake will execute the rule by invoking `cwltool`, which has to be available via your `$PATH` variable, and can be, e.g., installed via `conda` or `pip`.
When using in combination with :ref:`--software-deployment-method apptainer <apptainer>` (``--sdm`` for short), Snakemake will instruct `cwltool` to execute the command via Singularity in user space.
Otherwise, `cwltool` will in most cases use a Docker container, which requires Docker to be set up properly.

The advantage is that predefined tools available via any `repository of CWL tool definitions <https://www.commonwl.org/#Repositories_of_CWL_Tools_and_Workflows>`__ can be used in any supporting workflow management system.
In contrast to a :ref:`Snakemake wrapper <snakefiles-wrappers>`, CWL tool definitions are in general not suited to alter the behavior of a tool, e.g., by normalizing output names or special input handling.
As you can see in comparison to the analog :ref:`wrapper declaration <snakefiles-wrappers>` above, the rule becomes slightly more verbose, because input, output, and params have to be dispatched to the specific expectations of the CWL tool definition.



.. _snakefiles-modules-pathvars:

Pathvars
~~~~~~~~

It is possible to define :ref:`pathvars <snakefiles-pathvars>` on a per-module base as follows:

.. code-block:: python

    module other_workflow:
        snakefile:
            # here, plain paths, URLs and the special markers for code hosting providers (see below) are possible.
            "other_workflow/Snakefile"
        pathvars:
            results="results/other_workflow"

All rules in the module that make use of the defined pathvars will use whatever new values are defined there.
The given values will override eventual pathvar definitions inside of the module.

If further a config is passed to the module, any pathvar definitions in the config take precedence over pathvar definitions in the module definition, e.g.:

.. code-block:: python

    module other_workflow:
        snakefile:
            # here, plain paths, URLs and the special markers for code hosting providers (see below) are possible.
            "other_workflow/Snakefile"
        pathvars:
            results="results/other_workflow"
        config: config["other-workflow"]

Here, if ``config["other-workflow"]`` contains a ``pathvars`` section, those definitions will extend (and overwrite if also containing ``"results"``) the ``results`` pathvar defined in the module statement.
Concretely, consider the following two cases.
First, if the config contains the following:

.. code-block:: yaml

    pathvars:
        resources: "custom/path/to/resources"

Then inside of the module the two considered pathvars will be ``results="results/other_workflow"`` and ``resources="custom/path/to/resources"``.
If instead the config contains:
.. code-block:: yaml

    pathvars:
        results: "custom/results"
        resources: "custom/path/to/resources"

Then inside of the module the two considered pathvars will be ``results="custom/results"`` and ``resources="custom/path/to/resources"``.

Note that defining pathvars in the config should be considered a rare, discouraged and advanced use case, since the users has to know about the internal pathvar expectations of the module.
Workflow authors can explicitly forbid the modification of particular pathvars via :ref:`config file schemas and validation <snakefiles_config_validation>`.

.. _snakefiles-dynamic-modules:

Dynamic Modules
~~~~~~~~~~~~~~~

With Snakemake 9.0 and later, it is possible to load modules dynamically by providing the ``name`` keyword inside the module definition.
For example, by reading the module name from a config file or by iterating over several modules in a loop.
For this, the module name is not specified directly after the ``module`` keyword, but by specifying the ``name`` parameter.


.. code-block:: python

    for module_name in ['module1', 'module2']:
        module:
            name: module_name
            snakefile: f"{module_name}/Snakefile"
            config: config[module_name]

        use rule * from module_name as module_name*

.. note::
    It is not allowed to specify the module name both after the ``module`` keyword and inside the module definition after the ``name`` parameter.

In the ``use rule`` statement, it is first checked if the module name (here, ``'module_name'``) corresponds to a loaded module. If yes, the rules are imported from the loaded module and an arbitrary alias can be provided after the ``as`` keyword.

If ``module_name`` was not registered as a module (as in the example above), the module name is resolved dynamically by searching the name in the current python variable scope. In the example, it resolves to ``'module1'`` and ``'module2'``.
Note that this means that if ``use rule`` is used with the optional ``as`` keyword inside the loop, the alias after ``as`` must be specified using a variable to ensure a one-to-one mapping between module names and their aliases. This can either be the same name variable (as in the above example) or a second variable (as in the example below).

In particular, it is not possible to modify the alias name in the ``use rule`` statement (e.g., writing directly ``use rule * from module as module_*`` is not allowed for dynamic modules).

.. code-block:: python

    for module_name, alias in zip(['module1', 'module2'], ['module1_', 'module2_']):
        module:
            name: module_name
            snakefile: f"{module_name}/Snakefile"
            config: config[module_name]

        use rule * from module_name as alias*

..  _snakefiles-meta-wrappers:

Meta-Wrappers
~~~~~~~~~~~~~

Snakemake wrappers offer a simple way to include commonly used tools in Snakemake workflows.
In addition the `Snakemake Wrapper Repository`_ offers so-called meta-wrappers, which are combinations of wrappers, meant to perform common tasks.
Both wrappers and meta-wrappers are continuously tested.
The module statement also allows to easily use meta-wrappers, for example:

.. code-block:: python

    from snakemake.utils import min_version
    min_version("6.0")

    configfile: "config.yaml"


    module bwa_mapping:
        meta_wrapper: "0.72.0/meta/bio/bwa_mapping"


    use rule * from bwa_mapping


    def get_input(wildcards):
        return config["samples"][wildcards.sample]


    use rule bwa_mem from bwa_mapping with:
        input:
            get_input


First, we define the meta-wrapper as a module.
Next, we declare all rules from the module to be used.
And finally, we overwrite the input directive of the rule ``bwa_mem`` such that the raw data is taken from the place where our workflow configures it via it's config file.


.. _snakefile-code-hosting-providers:

Code hosting providers
~~~~~~~~~~~~~~~~~~~~~~

To obtain the correct URL to an external source code resource (e.g. a snakefile, see :ref:`snakefiles-modules`), Snakemake provides markers for code hosting providers.
Currently, Github

.. code-block:: python

    github("owner/repo", path="workflow/Snakefile", tag="v1.0.0")


and Gitlab are supported:

.. code-block:: python

    gitlab("owner/repo", path="workflow/Snakefile", tag="v1.0.0")

For the latter, it is also possible to specify an alternative host, e.g.

.. code-block:: python

    gitlab("owner/repo", path="workflow/Snakefile", tag="v1.0.0", host="somecustomgitlab.org")


While specifying a tag is highly encouraged, it is alternatively possible to specify a `commit` or a `branch` via respective keyword arguments.
Note that only when specifying a tag or a commit, Snakemake is able to persistently cache the source, thereby avoiding to repeatedly query it in case of multiple executions.


Private repositories
~~~~~~~~~~~~~~~~~~~~

To access source code resources located in private repositories you can define an
access token in the ``GITHUB_TOKEN`` and/or ``GITLAB_TOKEN`` environment variables.
