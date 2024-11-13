.. _tutorial-advanced:

Advanced: Decorating the example workflow
-----------------------------------------

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

Now that the basic concepts of Snakemake have been illustrated, we can introduce some advanced functionality.

Step 1: Specifying the number of used threads
:::::::::::::::::::::::::::::::::::::::::::::

For some tools, it is advisable to use more than one thread in order to speed up the computation.
**Snakemake can be made aware of the threads a rule needs** with the ``threads`` directive.
In our example workflow, it makes sense to use multiple threads for the rule ``bwa_map``:

.. code:: python

    rule bwa_map:
        input:
            "data/genome.fa",
            "data/samples/{sample}.fastq"
        output:
            "mapped_reads/{sample}.bam"
        threads: 8
        shell:
            "bwa mem -t {threads} {input} | samtools view -Sb - > {output}"

The number of threads can be propagated to the shell command with the familiar braces notation (i.e. ``{threads}``).
If no ``threads`` directive is given, a rule is assumed to need 1 thread.

When a workflow is executed, **the number of threads the jobs need is considered by the Snakemake scheduler**.
In particular, the scheduler ensures that the sum of the threads of all jobs running at the same time does not exceed a given number of available CPU cores.
This number is given with the ``--cores`` command line argument, which is mandatory for ``snakemake`` calls that actually run the workflow.
For example

.. code:: console

    $ snakemake --cores 10


.. note::

  Apart from the very common thread resource, Snakemake provides a ``resources`` directive that can be used to **specify arbitrary resources**, e.g., memory usage or auxiliary computing devices like GPUs.
  Similar to threads, these can be considered by the scheduler when an available amount of that resource is given with the command line argument ``--resources`` (see :ref:`snakefiles-resources`).

would execute the workflow with 10 cores.
Since the rule ``bwa_map`` needs 8 threads, only one job of the rule can run at a time, and the Snakemake scheduler will try to saturate the remaining cores with other jobs like, e.g., ``samtools_sort``.
The threads directive in a rule is interpreted as a maximum: when **less cores than threads** are provided, the number of threads a rule uses will be **reduced to the number of given cores**.

If ``--cores all`` is given, all available cores are used.

Exercise
........

* With the flag ``--forceall`` you can enforce a complete re-execution of the workflow. Combine this flag with different values for ``--cores`` and examine how the scheduler selects jobs to run in parallel.

Step 2: Config files
::::::::::::::::::::

So far, we specified which samples to consider by providing a Python list in the Snakefile.
However, often you want your workflow to be customizable, so that it can easily be adapted to new data.
For this purpose, Snakemake provides a `config file mechanism <https://snakemake.readthedocs.io/en/latest/snakefiles/configuration.html>`_.
Config files can be written in JSON_ or YAML_, and are used with the ``configfile`` directive.
In our example workflow, we add the line

.. code:: python

    configfile: "config.yaml"

to the top of the Snakefile.
Snakemake will load the config file and store its contents into a globally available `dictionary <https://docs.python.org/3/tutorial/datastructures.html#dictionaries>`_ named ``config``.
In our case, it makes sense to specify the samples in ``config.yaml`` as

.. code:: yaml

    samples:
        A: data/samples/A.fastq
        B: data/samples/B.fastq

Now, we can remove the statement defining ``SAMPLES`` from the Snakefile and change the rule ``bcftools_call`` to

.. code:: python

    rule bcftools_call:
        input:
            fa="data/genome.fa",
            bam=expand("sorted_reads/{sample}.bam", sample=config["samples"]),
            bai=expand("sorted_reads/{sample}.bam.bai", sample=config["samples"])
        output:
            "calls/all.vcf"
        shell:
            "bcftools mpileup -f {input.fa} {input.bam} | "
            "bcftools call -mv - > {output}"

.. _tutorial-input_functions:

Step 3: Input functions
:::::::::::::::::::::::

Since we have stored the path to the FASTQ files in the config file, we can also generalize the rule ``bwa_map`` to use these paths.
This case is different to the rule ``bcftools_call`` we modified above.
To understand this, it is important to know that Snakemake workflows are executed in three phases.

1. In the **initialization** phase, the files defining the workflow are parsed and all rules are instantiated.
2. In the **DAG** phase, the directed acyclic dependency graph of all jobs is built by filling wildcards and matching input files to output files.
3. In the **scheduling** phase, the DAG of jobs is executed, with jobs started according to the available resources.

The expand functions in the list of input files of the rule ``bcftools_call`` are executed during the initialization phase.
In this phase, we don't know about jobs, wildcard values and rule dependencies.
Hence, we cannot determine the FASTQ paths for rule ``bwa_map`` from the config file in this phase, because we don't even know which jobs will be generated from that rule.
Instead, we need to defer the determination of input files to the DAG phase.
This can be achieved by specifying an **input function** instead of a string as inside of the input directive.
For the rule ``bwa_map`` this works as follows:

.. code:: python

    def get_bwa_map_input_fastqs(wildcards):
        return config["samples"][wildcards.sample]
    
    rule bwa_map:
        input:
            "data/genome.fa",
            get_bwa_map_input_fastqs
        output:
            "mapped_reads/{sample}.bam"
        threads: 8
        shell:
            "bwa mem -t {threads} {input} | samtools view -Sb - > {output}"

Any normal function would work as well.
Input functions take as **single argument** a ``wildcards`` object, that allows to access the wildcards values via attributes (here ``wildcards.sample``).
They have to **return a string or a list of strings**, that are interpreted as paths to input files (here, we return the path that is stored for the sample in the config file).
Input functions are evaluated once the wildcard values of a job are determined.


Exercise
........

* In the ``data/samples`` folder, there is an additional sample ``C.fastq``. Add that sample to the config file and see how Snakemake wants to recompute the part of the workflow belonging to the new sample, when invoking with ``snakemake -n --forcerun bcftools_call``.

Step 4: Rule parameters
:::::::::::::::::::::::

Sometimes, shell commands are not only composed of input and output files and some static flags.
In particular, it can happen that additional parameters need to be set depending on the wildcard values of the job.
For this, Snakemake allows to **define arbitrary parameters** for rules with the ``params`` directive.
In our workflow, it is reasonable to annotate aligned reads with so-called read groups, that contain metadata like the sample name.
We modify the rule ``bwa_map`` accordingly:

.. code:: python

    rule bwa_map:
        input:
            "data/genome.fa",
            get_bwa_map_input_fastqs
        output:
            "mapped_reads/{sample}.bam"
        params:
            rg=r"@RG\tID:{sample}\tSM:{sample}"
        threads: 8
        shell:
            "bwa mem -R '{params.rg}' -t {threads} {input} | samtools view -Sb - > {output}"

.. note::

  The ``params`` directive can also take functions like in Step 3 to defer
  initialization to the DAG phase. In contrast to input functions, these can
  optionally take additional arguments ``input``, ``output``, ``threads``, and ``resources``.

Similar to input and output files, ``params`` can be accessed from the shell command, the Python based ``run`` block, or the script directive (see :ref:`tutorial-script`).

Exercise
........

* Variant calling can consider a lot of parameters. A particularly important one is the prior mutation rate (1e-3 per default). It is set via the flag ``-P`` of the ``bcftools call`` command. Consider making this flag configurable via adding a new key to the config file and using the ``params`` directive in the rule ``bcftools_call`` to propagate it to the shell command.

Step 5: Logging
:::::::::::::::

When executing a large workflow, it is usually desirable to store the logging output of each job into a separate file, instead of just printing all logging output to the terminal---when multiple jobs are run in parallel, this would result in chaotic output.
For this purpose, Snakemake allows to **specify log files** for rules.
Log files are defined via the ``log`` directive and handled similarly to output files, but they are not subject of rule matching and are not cleaned up when a job fails.
We modify our rule ``bwa_map`` as follows:

.. code:: python

    rule bwa_map:
        input:
            "data/genome.fa",
            get_bwa_map_input_fastqs
        output:
            "mapped_reads/{sample}.bam"
        params:
            rg=r"@RG\tID:{sample}\tSM:{sample}"
        log:
            "logs/bwa_mem/{sample}.log"
        threads: 8
        shell:
            "(bwa mem -R '{params.rg}' -t {threads} {input} | "
            "samtools view -Sb - > {output}) 2> {log}"

.. note::

  It is best practice to store all log files in a subdirectory ``logs/``, prefixed by the rule or tool name.

The shell command is modified to `collect STDERR output <https://tldp.org/LDP/abs/html/io-redirection.html>`_ of both ``bwa`` and ``samtools`` and pipe it into the file referred to by ``{log}``.
Log files must contain exactly the same wildcards as the output files to avoid file name clashes between different jobs of the same rule.

Exercise
........

* Add a log directive to the ``bcftools_call`` rule as well.
* Time to re-run the whole workflow (remember the command line flags to force re-execution). See how log files are created for variant calling and read mapping.
* The ability to track the provenance of each generated result is an important step towards reproducible analyses. Apart from the ``report`` functionality discussed before, Snakemake can summarize various provenance information for all output files of the workflow. The flag ``--summary`` prints a table associating each output file with the rule used to generate it, the creation date and optionally the version of the tool used for creation is provided. Further, the table informs about updated input files and changes to the source code of the rule after creation of the output file. Invoke Snakemake with ``--summary`` to examine the information for our example.

.. _tutorial_temp-and-protected-files:

Step 6: Temporary and protected files
:::::::::::::::::::::::::::::::::::::

In our workflow, we create two BAM files for each sample, namely
the output of the rules ``bwa_map`` and ``samtools_sort``.
When not dealing with examples, the underlying data is usually huge.
Hence, the resulting BAM files need a lot of disk space and their creation takes some time.
To save disk space, you can **mark output files as temporary**.
Snakemake will delete the marked files for you, once all the consuming jobs (that need it as input) have been executed.
We use this mechanism for the output file of the rule ``bwa_map``:

.. code:: python

    rule bwa_map:
        input:
            "data/genome.fa",
            get_bwa_map_input_fastqs
        output:
            temp("mapped_reads/{sample}.bam")
        params:
            rg=r"@RG\tID:{sample}\tSM:{sample}"
        log:
            "logs/bwa_mem/{sample}.log"
        threads: 8
        shell:
            "(bwa mem -R '{params.rg}' -t {threads} {input} | "
            "samtools view -Sb - > {output}) 2> {log}"

This results in the deletion of the BAM file once the corresponding ``samtools_sort`` job has been executed.
Since the creation of BAM files via read mapping and sorting is computationally expensive, it is reasonable to **protect** the final BAM file **from accidental deletion or modification**.
We modify the rule ``samtools_sort`` to mark its output file as ``protected``:

.. code:: python

    rule samtools_sort:
        input:
            "mapped_reads/{sample}.bam"
        output:
            protected("sorted_reads/{sample}.bam")
        shell:
            "samtools sort -T sorted_reads/{wildcards.sample} "
            "-O bam {input} > {output}"

After successful execution of the job, Snakemake will write-protect the output file in the filesystem, so that it can't be overwritten or deleted by accident.

Exercise
........

* Re-execute the whole workflow and observe how Snakemake handles the temporary and protected files.
* Run Snakemake with the target ``mapped_reads/A.bam``. Although the file is marked as temporary, you will see that Snakemake does not delete it because it is specified as a target file.
* Try to re-execute the whole workflow again with the dry-run option. You will see that it fails (as intended) because Snakemake cannot overwrite the protected output files.

After having a look at the summary, please go on with the :ref:`"additional features" part of the tutorial <tutorial-additional_features>`.

Summary
:::::::

For this advanced part of the tutorial, we have now created a ``config.yaml`` configuration file:

.. code:: yaml

    samples:
        A: data/samples/A.fastq
        B: data/samples/B.fastq
    
    prior_mutation_rate: 0.001


With this, the final version of our workflow in the ``Snakefile`` looks like this:

.. code:: python

    configfile: "config.yaml"


    rule all:
        input:
            "plots/quals.svg"


    def get_bwa_map_input_fastqs(wildcards):
        return config["samples"][wildcards.sample]
    

    rule bwa_map:
        input:
            "data/genome.fa",
            get_bwa_map_input_fastqs
        output:
            temp("mapped_reads/{sample}.bam")
        params:
            rg=r"@RG\tID:{sample}\tSM:{sample}"
        log:
            "logs/bwa_mem/{sample}.log"
        threads: 8
        shell:
            "(bwa mem -R '{params.rg}' -t {threads} {input} | "
            "samtools view -Sb - > {output}) 2> {log}"


    rule samtools_sort:
        input:
            "mapped_reads/{sample}.bam"
        output:
            protected("sorted_reads/{sample}.bam")
        shell:
            "samtools sort -T sorted_reads/{wildcards.sample} "
            "-O bam {input} > {output}"


    rule samtools_index:
        input:
            "sorted_reads/{sample}.bam"
        output:
            "sorted_reads/{sample}.bam.bai"
        shell:
            "samtools index {input}"


    rule bcftools_call:
        input:
            fa="data/genome.fa",
            bam=expand("sorted_reads/{sample}.bam", sample=config["samples"]),
            bai=expand("sorted_reads/{sample}.bam.bai", sample=config["samples"])
        output:
            "calls/all.vcf"
        params:
            rate=config["prior_mutation_rate"]
        log:
            "logs/bcftools_call/all.log"
        shell:
            "(bcftools mpileup -f {input.fa} {input.bam} | "
            "bcftools call -mv -P {params.rate} - > {output}) 2> {log}"


    rule plot_quals:
        input:
            "calls/all.vcf"
        output:
            "plots/quals.svg"
        script:
            "scripts/plot-quals.py"

Now, please go on with the :ref:`"additional features" part of the tutorial <tutorial-additional_features>`.