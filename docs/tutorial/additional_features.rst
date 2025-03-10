.. _tutorial-additional_features:

Additional features
-------------------

.. _Snakemake: https://snakemake.readthedocs.io
.. _Snakemake homepage: https://snakemake.readthedocs.io
.. _GNU Make: https://www.gnu.org/software/make
.. _Python: https://www.python.org
.. _BWA: http://bio-bwa.sourceforge.net
.. _SAMtools: https://www.htslib.org
.. _BCFtools: https://www.htslib.org
.. _Pandas: https://pandas.pydata.org
.. _Miniconda: https://conda.pydata.org/miniconda.html
.. _Conda: https://conda.pydata.org
.. _Bash: https://www.tldp.org/LDP/Bash-Beginners-Guide/html
.. _Atom: https://atom.io
.. _Anaconda: https://anaconda.org
.. _Graphviz: https://www.graphviz.org
.. _RestructuredText: https://docutils.sourceforge.io/docs/user/rst/quickstart.html
.. _data URI: https://developer.mozilla.org/en-US/docs/Web/HTTP/data_URIs
.. _JSON: https://json.org
.. _YAML: https://yaml.org
.. _DRMAA: https://www.drmaa.org
.. _rpy2: https://rpy2.github.io
.. _R: https://www.r-project.org
.. _Rscript: https://stat.ethz.ch/R-manual/R-devel/library/utils/html/Rscript.html
.. _PyYAML: https://pyyaml.org
.. _Docutils: https://docutils.sourceforge.io
.. _Bioconda: https://bioconda.github.io
.. _Vagrant: https://www.vagrantup.com
.. _Vagrant Documentation: https://docs.vagrantup.com
.. _Blogpost: https://blog.osteel.me/posts/2015/01/25/how-to-use-vagrant-on-windows.html
.. _slides: https://slides.com/johanneskoester/deck-1

In the following, we introduce some features that are beyond the scope of above example workflow.
For details and even more features, see :ref:`user_manual-writing_snakefiles`, :ref:`project_info-faq` and the command line help (``snakemake --help``).


Benchmarking
::::::::::::

With the ``benchmark`` directive, Snakemake can be instructed to **measure the wall clock time of a job**.
We activate benchmarking for the rule ``bwa_map``:

.. code:: python

    rule bwa_map:
        input:
            "data/genome.fa",
            lambda wildcards: config["samples"][wildcards.sample]
        output:
            temp("mapped_reads/{sample}.bam")
        params:
            rg="@RG\tID:{sample}\tSM:{sample}"
        log:
            "logs/bwa_mem/{sample}.log"
        benchmark:
            "benchmarks/{sample}.bwa.benchmark.txt"
        threads: 8
        shell:
            "(bwa mem -R '{params.rg}' -t {threads} {input} | "
            "samtools view -Sb - > {output}) 2> {log}"

The ``benchmark`` directive takes a string that points to the file where benchmarking results shall be stored.
Similar to output files, the path can contain wildcards (it must be the same wildcards as in the output files).
When a job derived from the rule is executed, Snakemake will measure the wall clock time and memory usage (in MiB) and store it in the file in tab-delimited format.
It is possible to repeat a benchmark multiple times in order to get a sense for the variability of the measurements.
This can be done by annotating the benchmark file, e.g., with ``repeat("benchmarks/{sample}.bwa.benchmark.txt", 3)`` Snakemake can be told to run the job three times.
The repeated measurements occur as subsequent lines in the tab-delimited benchmark file.

Modularization
::::::::::::::

In order to reuse building blocks or simply to structure large workflows, it is sometimes reasonable to **split a workflow into multiple Snakefiles**.
For this, Snakemake provides the ``include`` directive to include another Snakefile into the current one, e.g.:

.. code:: python

    include: "path/to/other.smk"

As can be seen, the default file extensions for snakefiles other than the main snakefile is ``.smk``.
Alternatively, Snakemake allows to **define external workflows as modules**.
A module refers to a working directory with a complete Snakemake workflow.
Output files of that sub-workflow can be used in the current Snakefile.
When executing, Snakemake ensures that the output files of the sub-workflow are up-to-date before executing the current workflow.
This mechanism is particularly useful when you want to extend a previous analysis without modifying it.
For details about modules, see the :ref:`documentation <snakefiles-modules>`.


Exercise
........

* Put the read mapping related rules into a separate Snakefile and use the ``include`` directive to make them available in our example workflow again.


.. _tutorial-conda:

Automatic deployment of software dependencies
:::::::::::::::::::::::::::::::::::::::::::::

In order to get a fully reproducible data analysis, it is not sufficient to
be able to execute each step and document all used parameters.
The used software tools and libraries have to be documented as well.
In this tutorial, you have already seen how Conda_ can be used to specify an
isolated software environment for a whole workflow. With Snakemake, you can
go one step further and specify Conda environments per rule.
This way, you can even make use of conflicting software versions (e.g. combine
Python 2 with Python 3).

In our example, instead of using an external environment we can specify
environments per rule, e.g.:

.. code:: python

  rule samtools_index:
    input:
        "sorted_reads/{sample}.bam"
    output:
        "sorted_reads/{sample}.bam.bai"
    conda:
        "envs/samtools.yaml"
    shell:
        "samtools index {input}"

with ``envs/samtools.yaml`` defined as

.. code:: yaml

  channels:
    - bioconda
    - conda-forge
  dependencies:
    - samtools =1.9

.. note::

  The conda directive does not work in combination with ``run`` blocks, because
  they have to share their Python environment with the surrounding snakefile.

When Snakemake is executed with

.. code:: console

  snakemake --software-deployment-method conda --cores 1

or the short form

.. code:: console

    snakemake --sdm conda -c 1

it will automatically create required environments and
activate them before a job is executed.
It is best practice to specify at least the `major and minor version <https://semver.org/>`_ of any packages
in the environment definition. Specifying environments per rule in this way has two
advantages.
First, the workflow definition also documents all used software versions.
Second, a workflow can be re-executed (without admin rights)
on a vanilla system, without installing any
prerequisites apart from Snakemake and Miniconda_.


Tool wrappers
:::::::::::::

In order to simplify the utilization of popular tools, Snakemake provides a
repository of so-called wrappers
(the `Snakemake wrapper repository <https://snakemake-wrappers.readthedocs.io>`_).
A wrapper is a short script that wraps (typically)
a command line application and makes it directly addressable from within Snakemake.
For this, Snakemake provides the ``wrapper`` directive that can be used instead of
``shell``, ``script``, or ``run``.
For example, the rule ``bwa_map`` could alternatively look like this:

.. code:: python

  rule bwa_mem:
    input:
        ref="data/genome.fa",
        sample=lambda wildcards: config["samples"][wildcards.sample]
    output:
        temp("mapped_reads/{sample}.bam")
    log:
        "logs/bwa_mem/{sample}.log"
    params:
        "-R '@RG\tID:{sample}\tSM:{sample}'"
    threads: 8
    wrapper:
        "0.15.3/bio/bwa/mem"

.. note::

  Updates to the Snakemake wrapper repository are automatically tested via
  `continuous integration <https://en.wikipedia.org/wiki/Continuous_integration>`_.

The wrapper directive expects a (partial) URL that points to a wrapper in the repository.
These can be looked up in the corresponding `database <https://snakemake-wrappers.readthedocs.io>`_.
The first part of the URL is a Git version tag. Upon invocation, Snakemake
will automatically download the requested version of the wrapper.
Furthermore, in combination with ``--software-deployment-method conda`` (see :ref:`tutorial-conda`),
the required software will be automatically deployed before execution.

.. _tutorial-cloud:

Cluster or cloud execution
::::::::::::::::::::::::::

Executing jobs on a cluster or in the cloud is supported by so-called executor plugins, which are distributed and documented via the `Snakemake plugin catalog <https://snakemake.github.io/snakemake-plugin-catalog>`_.

Constraining wildcards
::::::::::::::::::::::

Snakemake uses regular expressions to match output files to input files and determine dependencies between the jobs.
Sometimes it is useful to constrain the values a wildcard can have.
This can be achieved by adding a regular expression that describes the set of allowed wildcard values.
For example, the wildcard ``sample`` in the output file ``"sorted_reads/{sample}.bam"`` can be constrained to only allow alphanumeric sample names as ``"sorted_reads/{sample,[A-Za-z0-9]+}.bam"``.
Constraints may be defined per rule or globally using the ``wildcard_constraints`` keyword, as demonstrated in :ref:`snakefiles-wildcards`.
This mechanism helps to solve two kinds of ambiguity.

* It can help to avoid ambiguous rules, i.e. two or more rules that can be applied to generate the same output file. Other ways of handling ambiguous rules are described in the Section :ref:`snakefiles-ambiguous-rules`.
* It can help to guide the regular expression based matching so that wildcards are assigned to the right parts of a file name. Consider the output file ``{sample}.{group}.txt`` and assume that the target file is ``A.1.normal.txt``. It is not clear whether ``dataset="A.1"`` and ``group="normal"`` or ``dataset="A"`` and ``group="1.normal"`` is the right assignment. Here, constraining the dataset wildcard by ``{sample,[A-Z]+}.{group}`` solves the problem.

When dealing with ambiguous rules, it is best practice to first try to solve the ambiguity by using a proper file structure, for example, by separating the output files of different steps in different directories.

As follow-up to this tutorial, we recommend to have a look at the :ref:`interaction, visualization and reporting tutorial <interaction_visualization_reporting_tutorial>`, which focuses on Snakemake's ability to cover the last mile of data analysis, i.e., the generation of publication ready reports and interactive visualizations.