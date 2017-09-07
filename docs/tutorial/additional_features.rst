.. tutorial-additional_features:

Additional features
-------------------

.. _Snakemake: https://snakemake.bitbucket.io
.. _Snakemake homepage: https://snakemake.bitbucket.io
.. _GNU Make: https://www.gnu.org/software/make
.. _Python: http://www.python.org
.. _BWA: http://bio-bwa.sourceforge.net
.. _SAMtools: http://www.htslib.org
.. _BCFtools: http://www.htslib.org
.. _Pandas: http://pandas.pydata.org
.. _Miniconda: http://conda.pydata.org/miniconda.html
.. _Conda: http://conda.pydata.org
.. _Bash: http://www.tldp.org/LDP/Bash-Beginners-Guide/html
.. _Atom: https://atom.io
.. _Anaconda: https://anaconda.org
.. _Graphviz: http://www.graphviz.org
.. _RestructuredText: http://docutils.sourceforge.net/rst.html
.. _data URI: https://developer.mozilla.org/en-US/docs/Web/HTTP/data_URIs
.. _JSON: http://json.org
.. _YAML: http://yaml.org
.. _DRMAA: http://www.drmaa.org
.. _rpy2: http://rpy.sourceforge.net
.. _R: https://www.r-project.org
.. _Rscript: https://stat.ethz.ch/R-manual/R-devel/library/utils/html/Rscript.html
.. _PyYAML: http://pyyaml.org
.. _Docutils: http://docutils.sourceforge.net
.. _Bioconda: https://bioconda.github.io
.. _Vagrant: https://www.vagrantup.com
.. _Vagrant Documentation: https://docs.vagrantup.com
.. _Blogpost: http://blog.osteel.me/posts/2015/01/25/how-to-use-vagrant-on-windows.html
.. _slides: http://slides.com/johanneskoester/deck-1

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
With the command line flag ``--benchmark-repeats``, Snakemake can be instructed to perform repetitive measurements by executing benchmark jobs multiple times.
The repeated measurements occur as subsequent lines in the tab-delimited benchmark file.

We can include the benchmark results into our report:

.. code:: python

    rule report:
        input:
            T1="calls/all.vcf",
            T2=expand("benchmarks/{sample}.bwa.benchmark.txt", sample=config["samples"])
        output:
            "report.html"
        run:
            from snakemake.utils import report
            with open(input.T1) as vcf:
                n_calls = sum(1 for l in vcf if not l.startswith("#"))

            report("""
            An example variant calling workflow
            ===================================

            Reads were mapped to the Yeast
            reference genome and variants were called jointly with
            SAMtools/BCFtools.

            This resulted in {n_calls} variants (see Table T1_).
            Benchmark results for BWA can be found in the tables T2_.
            """, output[0], **input)

We use the ``expand`` function to collect the benchmark files for all samples.
Here, we directly provide names for the input files.
In particular, we can also name the whole list of benchmark files returned by the ``expand`` function as ``T2``.
When invoking the ``report`` function, we just unpack ``input`` into keyword arguments (resulting in ``T1`` and ``T2``).
In the text, we refer with ``T2_`` to the list of benchmark files.

Exercise
........

* Re-execute the workflow and benchmark ``bwa_map`` with 3 repeats. Open the report and see how the list of benchmark files is presented in the HTML report.

Modularization
::::::::::::::

In order to re-use building blocks or simply to structure large workflows, it is sometimes reasonable to **split a workflow into modules**.
For this, Snakemake provides the ``include`` directive to include another Snakefile into the current one, e.g.:

.. code:: python

    include: "path/to/other.snakefile"

Alternatively, Snakemake allows to **define sub-workflows**.
A sub-workflow refers to a working directory with a complete Snakemake workflow.
Output files of that sub-workflow can be used in the current Snakefile.
When executing, Snakemake ensures that the output files of the sub-workflow are up-to-date before executing the current workflow.
This mechanism is particularly useful when you want to extend a previous analysis without modifying it.
For details about sub-workflows, see the :ref:`documentation <snakefiles-sub_workflows>`.


Exercise
........

* Put the read mapping related rules into a separate Snakefile and use the ``include`` directive to make them available in our example workflow again.


Using custom scripts
::::::::::::::::::::

Using the ``run`` directive as above is only reasonable for short Python scripts.
As soon as your script becomes larger, it is reasonable to separate it from the
workflow definition.
For this purpose, Snakemake offers the ``script`` directive.
Using this, the ``report`` rule from above could instead look like this:

.. code:: python

    rule report:
        input:
            T1="calls/all.vcf",
            T2=expand("benchmarks/{sample}.bwa.benchmark.txt", sample=config["samples"])
        output:
            "report.html"
        script:
            "scripts/report.py"

The actual Python code to generate the report is now hidden in the script ``scripts/report.py``.
Script paths are always relative to the referring Snakefile.
In the script, all properties of the rule like ``input``, ``output``, ``wildcards``,
``params``, ``threads`` etc. are available as attributes of a global ``snakemake`` object:

.. code:: python

    from snakemake.utils import report

    with open(snakemake.input.T1) as vcf:
        n_calls = sum(1 for l in vcf if not l.startswith("#"))

    report("""
    An example variant calling workflow
    ===================================

    Reads were mapped to the Yeast
    reference genome and variants were called jointly with
    SAMtools/BCFtools.

    This resulted in {n_calls} variants (see Table T1_).
    Benchmark results for BWA can be found in the tables T2_.
    """, snakemake.output[0], **snakemake.input)

.. sidebar:: Note

  It is best practice to use the script directive whenever a run block would have
  more than a few lines of code.

Although there are other strategies to invoke separate scripts from your workflow
(e.g., invoking them via shell commands), the benefit of this is obvious:
the script logic is separated from the workflow logic (and can be even shared between workflows),
but boilerplate code like the parsing of command line arguments in unnecessary.

Apart from Python scripts, it is also possible to use R scripts. In R scripts,
an S4 object named ``snakemake`` analog to the Python case above is available and
allows access to input and output files and other parameters. Here the syntax
follows that of S4 classes with attributes that are R lists, e.g. we can access
the first input file with ``snakemake@input[[1]]`` (note that the first file does
not have index 0 here, because R starts counting from 1). Named input and output
files can be accessed in the same way, by just providing the name instead of an
index, e.g. ``snakemake@input[["myfile"]]``.

For details and examples, see the :ref:`snakefiles-external_scripts` section in the Documentation.

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
  dependencies:
    - samtools =1.3

.. sidebar:: Note

  The conda directive does not work in combination with ``run`` blocks, because
  they have to share their Python environment with the surrounding snakefile.

When Snakemake is executed with

.. code:: console

  snakemake --use-conda

it will automatically create required environments and
activate them before a job is executed.
It is best practice to specify at least the `major and minor version <http://semver.org/>`_ of any packages
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

.. sidebar:: Note

  Updates to the Snakemake wrapper repository are automatically tested via
  `continuous integration <https://en.wikipedia.org/wiki/Continuous_integration>`_.

The wrapper directive expects a (partial) URL that points to a wrapper in the repository.
These can be looked up in the corresponding `database <https://snakemake-wrappers.readthedocs.io>`_.
The first part of the URL is a Git version tag. Upon invocation, Snakemake
will automatically download the requested version of the wrapper.
Furthermore, in combination with ``--use-conda`` (see :ref:`tutorial-conda`),
the required software will be automatically deployed before execution.

Cluster execution
:::::::::::::::::

By default, Snakemake executes jobs on the local machine it is invoked on.
Alternatively, it can execute jobs in **distributed environments, e.g., compute clusters or batch systems**.
If the nodes share a common file system, Snakemake supports three alternative execution modes.

In cluster environments, compute jobs are usually submitted as shell scripts via commands like ``qsub``.
Snakemake provides a **generic mode** to execute on such clusters.
By invoking Snakemake with

.. code:: console

    $ snakemake --cluster qsub --jobs 100

each job will be compiled into a shell script that is submitted with the given command (here ``qsub``).
The ``--jobs`` flag limits the number of concurrently submitted jobs to 100.
This basic mode assumes that the submission command returns immediately after submitting the job.
Some clusters allow to run the submission command in **synchronous mode**, such that it waits until the job has been executed.
In such cases, we can invoke e.g.

.. code:: console

    $ snakemake --cluster-sync "qsub -sync yes" --jobs 100

The specified submission command can also be **decorated with additional parameters taken from the submitted job**.
For example, the number of used threads can be accessed in braces similarly to the formatting of shell commands, e.g.

.. code:: console

    $ snakemake --cluster "qsub -pe threaded {threads}" --jobs 100

Alternatively, Snakemake can use the Distributed Resource Management Application API (DRMAA_).
This API provides a common interface to control various resource management systems.
The **DRMAA support** can be activated by invoking Snakemake as follows:

.. code:: console

    $ snakemake --drmaa --jobs 100

If available, **DRMAA is preferable over the generic cluster modes** because it provides better control and error handling.
To support additional cluster specific parametrization, a Snakefile can be complemented by a :ref:`snakefiles-cluster_configuration` file.


Constraining wildcards
::::::::::::::::::::::

Snakemake uses regular expressions to match output files to input files and determine dependencies between the jobs.
Sometimes it is useful to constrain the values a wildcard can have.
This can be achieved by adding a regular expression that describes the set of allowed wildcard values.
For example, the wildcard ``sample`` in the output file ``"sorted_reads/{sample}.bam"`` can be constrained to only allow alphanumeric sample names as ``"sorted_reads/{sample,[A-Za-z0-9]+}.bam"``.
Constrains may be defined per rule or globally using the ``wildcard_constraints`` keyword, as demonstrated in :ref:`snakefiles-wildcards`.
This mechanism helps to solve two kinds of ambiguity.

* It can help to avoid ambiguous rules, i.e. two or more rules that can be applied to generate the same output file. Other ways of handling ambiguous rules are described in the Section :ref:`snakefiles-ambiguous-rules`.
* It can help to guide the regular expression based matching so that wildcards are assigned to the right parts of a file name. Consider the output file ``{sample}.{group}.txt`` and assume that the target file is ``A.1.normal.txt``. It is not clear whether ``dataset="A.1"`` and ``group="normal"`` or ``dataset="A"`` and ``group="1.normal"`` is the right assignment. Here, constraining the dataset wildcard by ``{sample,[A-Z]+}.{group}`` solves the problem.

When dealing with ambiguous rules, it is best practice to first try to solve the ambiguity by using a proper file structure, for example, by separating the output files of different steps in different directories.
