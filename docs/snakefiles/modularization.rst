.. snakefiles-modularization:

.. _Snakemake Wrapper Repository: https://snakemake-wrappers.readthedocs.io

==============
Modularization
==============

Modularization in Snakemake comes at different levels.

1. The most fine-grained level are wrappers. They are available and can be published at the `Snakemake Wrapper Repository`_. These wrappers can then be composed and customized according to your needs, by copying skeleton rules into your workflow. In combination with conda integration, wrappers also automatically deploy the needed software dependencies into isolated environments.
2. For larger, reusable parts that shall be integrated into a common workflow, it is recommended to write small Snakefiles and include them into a master Snakefile via the include statement. In such a setup, all rules share a common config file.
3. The third level of separation are subworkflows. Importantly, these are rather meant as links between otherwise separate data analyses.


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
Thereby, 0.0.8 can be replaced with the git `version tag <https://github.com/snakemake/snakemake-wrappers/releases>`_ you want to use, or a `commit id <https://github.com/snakemake/snakemake-wrappers/commits>`_.
This ensures reproducibility since changes in the wrapper implementation won't be propagated automatically to your workflow.
Alternatively, e.g., for development, the wrapper directive can also point to full URLs, including URLs to local files with absolute paths ``file://`` or relative paths ``file:``.
Examples for each wrapper can be found in the READMEs located in the wrapper subdirectories at the `Snakemake Wrapper Repository`_.

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
