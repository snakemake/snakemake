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

.. _snakefiles-includes:

--------
Includes
--------

Another Snakefile with all its rules can be included into the current:

.. code-block:: python

    include: "path/to/other/snakefile"

The default target rule (often called the ``all``-rule), won't be affected by the include.
I.e. it will always be the first rule in your Snakefile, no matter how many includes you have above your first rule.
Includes are relative to the directory of the Snakefile in which they occur.
For example, if above Snakefile resides in the directory ``my/dir``, then Snakemake will search for the include at ``my/dir/path/to/other/snakefile``, regardless of the working directory.


.. _snakefiles-modules:

-------
Modules
-------

With Snakemake 6.0 and later, it is possible to define external workflows as modules, from which rules can be used by explicitly "importing" them.

.. code-block:: python

    from snakemake.utils import min_version
    min_version("6.0")

    module other_workflow:
        snakefile:
            # here, plain paths, URLs and the special markers for code hosting providers (see below) are possible.
            "other_workflow/Snakefile"
    
    use rule * from other_workflow exclude ruleC as other_*

The ``module other_workflow:`` statement registers the external workflow as a module, by defining the path to the main snakefile of ``other_workflow``.
Here, plain paths, HTTP/HTTPS URLs and special markers for code hosting providers like Github or Gitlab are possible (see :ref:`snakefile-code-hosting-providers`).
The second statement, ``use rule * from other_workflow exclude ruleC as other_*``, declares all rules of that module to be used in the current one, except for ruleC.
Thereby, the ``as other_*`` at the end renames all those rules with a common prefix.
This can be handy to avoid rule name conflicts (note that rules from modules can otherwise overwrite rules from your current workflow or other modules).

.. note::

    The imported module cannot be named as `workflow`, which is a reserved name.

The module is evaluated in a separate namespace, and only the selected rules are added to the current workflow.
Non-rule Python statements inside the module are also evaluated in that separate namespace.
They are available in the module-defining workflow under the name of the module (e.g. here ``other_workflow.myfunction()`` would call the function ``myfunction`` that has been defined in the model, e.g. in ``other_workflow/Snakefile``).
Also note that this means that any Python variables and functions available in the module-defining namespace will **not** be visible from inside the module.
However, it is possible to pass information to the module using the ``config`` mechanism described in the following.

It is possible to overwrite the global config dictionary for the module, which is usually filled by the ``configfile`` statement (see :ref:`snakefiles_standard_configuration`):

.. code-block:: python

    from snakemake.utils import min_version
    min_version("6.0")

    configfile: "config/config.yaml"

    module other_workflow:
        # here, plain paths, URLs and the special markers for code hosting providers (see below) are possible.
        snakefile: "other_workflow/Snakefile"
        config: config["other-workflow"]
    
    use rule * from other_workflow as other_*

In this case, any ``configfile`` statements inside the module are ignored.
In addition, it is possible to skip any :ref:`validation <snakefiles_config_validation>` statements in the module, by specifying ``skip_validation: True`` in the module statement.
Moreover, one can automatically move all relative input and output files of a module into a dedicated folder by specifying ``prefix: "foo"`` in the module definition, e.g. any output file ``path/to/output.txt`` in the module would be stored under ``foo/path/to/output.txt`` instead.
This becomes particularly useful when combining multiple modules, see :ref:`use_with_modules`.
However, if you have some input files that come from outside the workflow, you can use the ``local`` flag so that their path is not modified (see :ref:`snakefiles-storage-local-files`)..

Instead of using all rules, it is possible to import specific rules.
Specific rules may even be modified before using them, via a final ``with:`` followed by a block that lists items to overwrite.
This modification can be performed after a general import, and will overwrite any unmodified import of the same rule.

.. code-block:: python

    from snakemake.utils import min_version
    min_version("6.0")

    module other_workflow:
        # here, plain paths, URLs and the special markers for code hosting providers (see below) are possible.
        snakefile: "other_workflow/Snakefile"
        config: config["other-workflow"]

    use rule * from other_workflow as other_*

    use rule some_task from other_workflow as other_some_task with:
        output:
            "results/some-result.txt"

By such a modifying use statement, any properties of the rule (``input``, ``output``, ``log``, ``params``, ``benchmark``, ``threads``, ``resources``, etc.) can be overwritten, except the actual execution step (``shell``, ``notebook``, ``script``, ``cwl``, or ``run``).

Note that the second use statement has to use the original rule name, not the one that has been prefixed with ``other_`` via the first use statement (there is no rule ``other_some_task`` in the module ``other_workflow``).
In order to overwrite the rule ``some_task`` that has been imported with the first ``use rule`` statement, it is crucial to ensure that the rule is used with the same name in the second statement, by adding an equivalent ``as`` clause (here ``other_some_task``).
Otherwise, you will have two versions of the same rule, which might be unintended (a common symptom of such unintended repeated uses would be ambiguous rule exceptions thrown by Snakemake).

Of course, it is possible to combine the use of rules from multiple modules (see :ref:`use_with_modules`), and via modifying statements they can be rewired and reconfigured in an arbitrary way.

..  _snakefiles-meta-wrappers:

~~~~~~~~~~~~~
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

----------------------
Code hosting providers
----------------------

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

~~~~~~~~~~~~~~~~~~~~
Private repositories
~~~~~~~~~~~~~~~~~~~~

To access source code resources located in private repositories you can define an
access token in the ``GITHUB_TOKEN`` and/or ``GITLAB_TOKEN`` environment variables.
