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

The `wrapper` directive allows to have reusable wrapper scripts around e.g. command line tools.

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


Source files can also be provided as plain HTTP/HTTPS URLs.
In that case, they are treated as generic remote source files.

As a convenience syntax, hosted source files can also be written with an explicit provider prefix.
Use ``gh:owner/repo@ref:path/to/Snakefile`` for GitHub, optionally inserting a custom host as ``gh:github.example.org:owner/repo@ref:path/to/Snakefile``.
Use ``gl:group/project@ref:path/to/Snakefile`` for GitLab, optionally inserting a custom host as ``gl:gitlab.example.org:group/project@ref:path/to/Snakefile``.
If the trailing path is omitted, Snakemake assumes ``workflow/Snakefile``.
When refs contain slashes, prefer this shorthand because it avoids ambiguity between ref and file path.
For example:

.. code-block:: python

    "gh:snakemake-workflows/dna-seq-gatk-variant-calling@v2.0.1:workflow/Snakefile"
    "gh:snakemake-workflows/dna-seq-gatk-variant-calling@v2.0.1"
    "gh:github.example.org:owner/repo@main:workflow/Snakefile"
    "gl:owner/repo@main:workflow/Snakefile"
    "gl:gitlab.cern.ch:group/project@main:workflow/Snakefile"

While specifying a tag is highly encouraged, it is alternatively possible to specify a `commit` or a `branch` via respective keyword arguments.
Note that only when specifying a tag or a commit, Snakemake is able to persistently cache the source, thereby avoiding to repeatedly query it in case of multiple executions.


Private repositories
~~~~~~~~~~~~~~~~~~~~

To access source code resources located in private repositories you can define an
access token in the ``GITHUB_TOKEN`` and/or ``GITLAB_TOKEN`` environment variables.
