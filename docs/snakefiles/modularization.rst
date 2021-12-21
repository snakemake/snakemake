.. snakefiles-modularization:

.. _Snakemake Wrapper Repository: https://snakemake-wrappers.readthedocs.io

==============
Modularization
==============

Modularization in Snakemake comes at four different levels.

1. The most fine-grained level are wrappers. They are available and can be published at the `Snakemake Wrapper Repository`_. These wrappers can then be composed and customized according to your needs, by copying skeleton rules into your workflow. In combination with conda integration, wrappers also automatically deploy the needed software dependencies into isolated environments.
2. For larger, reusable parts that shall be integrated into a common workflow, it is recommended to write small Snakefiles and include them into a main Snakefile via the include statement. In such a setup, all rules share a common config file.
3. The third level is provided via the :ref:`module statement <snakefiles-modules>`, which enables arbitrary combination and reuse of rules.
4. Finally, Snakemake provides a syntax for defining :ref:`subworkflows <snakefiles-sub_workflows>`, which is however deprecated in favor of the module statement.


.. _snakefiles-wrappers:

--------
Wrappers
--------

The wrapper directive allows to have re-usable wrapper scripts around e.g. command line tools.
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

With Snakemake 4.8.0, it is possible to refer to `CWL <https://www.commonwl.org/>`_ tool definitions in rules instead of specifying a wrapper or a plain shell command.
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
When using in combination with :ref:`--use-singularity <singularity>`, Snakemake will instruct `cwltool` to execute the command via Singularity in user space.
Otherwise, `cwltool` will in most cases use a Docker container, which requires Docker to be set up properly.

The advantage is that predefined tools available via any `repository of CWL tool definitions <https://www.commonwl.org/#Repositories_of_CWL_Tools_and_Workflows>`_ can be used in any supporting workflow management system.
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
    
    use rule * from other_workflow as other_*

The first statement registers the external workflow as a module, by defining the path to the main snakefile.
Here, plain paths, HTTP/HTTPS URLs and special markers for code hosting providers like Github or Gitlab are possible (see :ref:`snakefile-code-hosting-providers`).
The second statement declares all rules of that module to be used in the current one.
Thereby, the ``as other_*`` at the end renames all those rule with a common prefix.
This can be handy to avoid rule name conflicts (note that rules from modules can otherwise overwrite rules from your current workflow or other modules).
The module is evaluated in a separate namespace, and only the selected rules are added to the current workflow.
Non-rule Python statements inside the module are also evaluated in that separate namespace.
They are available in the module-defining workflow under the name of the module (e.g. here ``other_workflow.myfunction()`` would call the function ``myfunction`` that has been define in the model, e.g. in ``other_workflow/Snakefile``).

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
In addition, it is possible to skip any :ref:`validation <snakefiles_config_validation>` statements in the module, by specifying ``skip_validation: True`` in the module statment.
Moreover, one can automatically move all relative input and output files of a module into a dedicated folder: by specifying ``prefix: "foo"`` in the module definition, e.g. any output file ``path/to/output.txt`` in the module would be stored under ``foo/path/to/output.txt`` instead.
This becomes particularly usefull when combining multiple modules, see :ref:`use_with_modules`.

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
Both wrappers and meta-wrappers are continously tested.
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

.. _snakefiles-sub_workflows:

-------------
Sub-Workflows
-------------

In addition to including rules of another workflow, Snakemake allows to depend on the output of other workflows as sub-workflows.
A sub-workflow is executed independently before the current workflow is executed.
Thereby, Snakemake ensures that all files the current workflow depends on are created or updated if necessary.
This allows to create links between otherwise separate data analyses.

.. code-block:: python

    subworkflow otherworkflow:
        workdir:
            "../path/to/otherworkflow"
        snakefile:
            "../path/to/otherworkflow/Snakefile"
        configfile:
            "path/to/custom_configfile.yaml"

    rule a:
        input:
            otherworkflow("test.txt")
        output: ...
        shell:  ...

Here, the subworkflow is named "otherworkflow" and it is located in the working directory ``../path/to/otherworkflow``.
The snakefile is in the same directory and called ``Snakefile``.
If ``snakefile`` is not defined for the subworkflow, it is assumed be located in the workdir location and called ``Snakefile``, hence, above we could have left the ``snakefile`` keyword out as well.
If ``workdir`` is not specified, it is assumed to be the same as the current one.
The (optional) definition of a ``configfile`` allows to parameterize the subworkflow as needed.
Files that are output from the subworkflow that we depend on are marked with the ``otherworkflow`` function (see the input of rule a).
This function automatically determines the absolute path to the file (here ``../path/to/otherworkflow/test.txt``).

When executing, snakemake first tries to create (or update, if necessary) ``test.txt`` (and all other possibly mentioned dependencies) by executing the subworkflow.
Then the current workflow is executed.
This can also happen recursively, since the subworkflow may have its own subworkflows as well.


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
