.. tutorial-advanced:

Advanced: Decorating the example workflow
-----------------------------------------

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

Now that the basic concepts of Snakemake have been illustrated, we can introduce advanced topics.

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
In particular, the scheduler ensures that the sum of the threads of all running jobs does not exceed a given number of available CPU cores.
This number can be given with the ``--cores`` command line argument (per default, Snakemake uses only 1 CPU core).
For example

.. code:: console

    $ snakemake --cores 10


.. sidebar:: Note

  Apart from the very common thread resource, Snakemake provides a ``resources`` directive that can be used to **specify arbitrary resources**, e.g., memory usage or auxiliary computing devices like GPUs.
  Similar to threads, these can be considered by the scheduler when an available amount of that resource is given with the command line argument ``--resources`` (see :ref:`snakefiles-resources`).

would execute the workflow with 10 cores.
Since the rule ``bwa_map`` needs 8 threads, only one job of the rule can run at a time, and the Snakemake scheduler will try to saturate the remaining cores with other jobs like, e.g., ``samtools_sort``.
The threads directive in a rule is interpreted as a maximum: when **less cores than threads** are provided, the number of threads a rule uses will be **reduced to the number of given cores**.

Exercise
........

* With the flag ``--forceall`` you can enforce a complete re-execution of the workflow. Combine this flag with different values for ``--cores`` and examine how the scheduler selects jobs to run in parallel.

Step 2: Config files
::::::::::::::::::::

So far, we specified the samples to consider in a Python list within the Snakefile.
However, often you want your workflow to be customizable, so that it can be easily adapted to new data.
For this purpose, Snakemake provides a config file mechanism.
Config files can be written in JSON_ or YAML_, and loaded with the ``configfile`` directive.
In our example workflow, we add the line

.. code:: python

    configfile: "config.yaml"

to the top of the Snakefile.
Snakemake will load the config file and store its contents into a globally available dictionary named ``config``.
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
            "samtools mpileup -g -f {input.fa} {input.bam} | "
            "bcftools call -mv - > {output}"

.. _tutorial-input_functions:

Step 3: Input functions
:::::::::::::::::::::::

Since we have stored the path to the FASTQ files in the config file, we can also generalize the rule ``bwa_map`` to use these paths.
This case is different to the rule ``bcftools_call`` we modified above.
To understand this, it is important to know that Snakemake workflows are executed in three phases.

* In the **initialization** phase, the workflow is parsed and all rules are instantiated.
* In the **DAG** phase, the DAG of jobs is built by filling wildcards and matching input files to output files.
* In the **scheduling** phase, the DAG of jobs is executed.

The expand functions in the list of input files of the rule ``bcftools_call`` are executed during the initialization phase.
In this phase, we don't know about jobs, wildcard values and rule dependencies.
Hence, we cannot determine the FASTQ paths for rule ``bwa_map`` from the config file in this phase, because we don't even know which jobs will be generated from that rule.
Instead, we need to defer the determination of input files to the DAG phase.
This can be achieved by specifying an **input function** instead of a string as inside of the input directive.
For the rule ``bwa_map`` this works as follows:

.. code:: python

    rule bwa_map:
        input:
            "data/genome.fa",
            lambda wildcards: config["samples"][wildcards.sample]
        output:
            "mapped_reads/{sample}.bam"
        threads: 8
        shell:
            "bwa mem -t {threads} {input} | samtools view -Sb - > {output}"

.. sidebar:: Note

  Snakemake does not automatically rerun jobs when new input files are added as
  in the excercise below. However, you can get a list of output files that
  are affected by such changes with ``snakemake --list-input-changes``.
  To trigger a rerun, some bash magic helps:

  .. code:: console

    snakemake -n --forcerun $(snakemake --list-input-changes)

Here, we use an anonymous function, also called `lambda expression <https://docs.python.org/3/tutorial/controlflow.html#lambda-expressions>`_.
Any normal function would work as well.
Input functions take as **single argument** a ``wildcards`` object, that allows to access the wildcards values via attributes (here ``wildcards.sample``).
They have to **return a string or a list of strings**, that are interpreted as paths to input files (here, we return the path that is stored for the sample in the config file).
Input functions are evaluated once the wildcard values of a job are determined.


Exercise
........

* In the ``data/samples`` folder, there is an additional sample ``C.fastq``. Add that sample to the config file and see how Snakemake wants to recompute the part of the workflow belonging to the new sample, when invoking with ``snakemake -n --reason --forcerun bcftools_call``.

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
            lambda wildcards: config["samples"][wildcards.sample]
        output:
            "mapped_reads/{sample}.bam"
        params:
            rg=r"@RG\tID:{sample}\tSM:{sample}"
        threads: 8
        shell:
            "bwa mem -R '{params.rg}' -t {threads} {input} | samtools view -Sb - > {output}"

.. sidebar:: Note

  The ``params`` directive can also take functions like in Step 3 to defer
  initialization to the DAG phase. In contrast to input functions, these can
  optionally take additional arguments ``input``, ``output``, ``threads``, and ``resources``.

Similar to input and output files, ``params`` can be accessed from the shell command the Python based ``run`` block, or the script directive (see :ref:`tutorial-script`).

Exercise
........

* Variant calling can consider a lot of parameters. A particularly important one is the prior mutation rate (1e-3 per default). It is set via the flag ``-P`` of the ``bcftools call`` command. Consider making this flag configurable via adding a new key to the config file and using the ``params`` directive in the rule ``bcftools_call`` to propagate it to the shell command.

Step 5: Logging
:::::::::::::::

When executing a large workflow, it is usually desirable to store the output of each job persistently in files instead of just printing it to the terminal.
For this purpose, Snakemake allows to **specify log files** for rules.
Log files are defined via the ``log`` directive and handled similarly to output files, but they are not subject of rule matching and are not cleaned up when a job fails.
We modify our rule ``bwa_map`` as follows:

.. code:: python

    rule bwa_map:
        input:
            "data/genome.fa",
            lambda wildcards: config["samples"][wildcards.sample]
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

.. sidebar:: Note

  It is best practice to store all log files in a subdirectory ``logs/``, prefixed by the rule or tool name.

The shell command is modified to collect STDERR output of both ``bwa`` and ``samtools`` and pipe it into the file referred by ``{log}``.
Log files must contain exactly the same wildcards as the output files to avoid clashes.

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
Snakemake allows to **mark output files as temporary**, such that they are deleted once every consuming job has been executed, in order to save disk space.
We use this mechanism for the output file of the rule ``bwa_map``:

.. code:: python

    rule bwa_map:
        input:
            "data/genome.fa",
            lambda wildcards: config["samples"][wildcards.sample]
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
We modify the rule ``samtools_sort`` by marking its output file as ``protected``:

.. code:: python

    rule samtools_sort:
        input:
            "mapped_reads/{sample}.bam"
        output:
            protected("sorted_reads/{sample}.bam")
        shell:
            "samtools sort -T sorted_reads/{wildcards.sample} "
            "-O bam {input} > {output}"

After execution of the job, Snakemake will write-protect the output file in the filesystem, so that it can't be overwritten or deleted accidentally.

Exercise
........

* Re-execute the whole workflow and observe how Snakemake handles the temporary and protected files.
* Run Snakemake with the target ``mapped_reads/A.bam``. Although the file is marked as temporary, you will see that Snakemake does not delete it because it is specified as a target file.
* Try to re-execute the whole workflow again with the dry-run option. You will see that it fails (as intended) because Snakemake cannot overwrite the protected output files.

Summary
:::::::

The final version of our workflow looks like this:

.. code:: python

    configfile: "config.yaml"


    rule all:
        input:
            "plots/quals.svg"


    rule bwa_map:
        input:
            "data/genome.fa",
            lambda wildcards: config["samples"][wildcards.sample]
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
        shell:
            "samtools mpileup -g -f {input.fa} {input.bam} | "
            "bcftools call -mv - > {output}"


    rule plot_quals:
        input:
            "calls/all.vcf"
        output:
            "plots/quals.svg"
        script:
            "scripts/plot-quals.py"
