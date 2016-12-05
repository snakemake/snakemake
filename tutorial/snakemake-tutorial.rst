Snakemake Tutorial
==================

.. _Snakemake: http://snakemake.bitbucket.org
.. _Snakemake homepage: http://snakemake.bitbucket.org
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
.. _Documentation: https://bitbucket.org/snakemake/snakemake/wiki/Documentation
.. _JSON: http://json.org
.. _YAML: http://yaml.org
.. _DRMAA: http://www.drmaa.org
.. _FAQ: https://bitbucket.org/snakemake/snakemake/wiki/FAQ
.. _rpy2: http://rpy.sourceforge.net
.. _R: https://www.r-project.org
.. _Rscript: https://stat.ethz.ch/R-manual/R-devel/library/utils/html/Rscript.html
.. _cluster configuration: https://bitbucket.org/snakemake/snakemake/wiki/Documentation#markdown-header-cluster-configuration
.. _script section in the Documentation: https://bitbucket.org/snakemake/snakemake/wiki/Documentation#markdown-header-external-scripts
.. _PyYAML: http://pyyaml.org
.. _Docutils: http://docutils.sourceforge.net
.. _Bioconda: https://bioconda.github.io
.. _Vagrant: https://www.vagrantup.com
.. _Vagrant Documentation: https://docs.vagrantup.com
.. _Blogpost: http://blog.osteel.me/posts/2015/01/25/how-to-use-vagrant-on-windows.html
.. _slides: http://slides.com/johanneskoester/deck-1

This tutorial introduces the text-based workflow system Snakemake_.
Snakemake follows the `GNU Make`_ paradigm: workflows are defined in terms of rules that define how to create output files from input files.
Dependencies between the rules are determined automatically, creating a DAG (directed acyclic graph) of jobs that can be automatically parallelized.

Snakemake sets itself apart from existing text-based workflow systems in the following way.
Hooking into the Python interpreter, Snakemake offers a definition language that is an extension of Python_ with syntax to define rules and workflow specific properties.
This allows to combine the flexibility of a plain scripting language with a pythonic workflow definition.
The Python language is known to be concise yet readable and can appear almost like pseudo-code.
The syntactic extensions provided by Snakemake maintain this property for the definition of the workflow.
Further, Snakemakes scheduling algorithm can be constrained by priorities, provided cores and customizable resources and it provides a generic support for distributed computing (e.g., cluster or batch systems).
Hence, a Snakemake workflow scales without modification from single core workstations and multi-core servers to cluster or batch systems.

While the examples presented here come from Bioinformatics, Snakemake is considered a general-purpose workflow management system for any discipline.

**Note:** To get the most up-to-date version of this tutorial and further information on Snakemake, please visit the `Snakemake homepage`_.
Also have a look at the corresponding slides_.

.. contents::

------------

Setup
-----

To go through this tutorial, you need the following software installed:

* Python_ ≥3.3
* Snakemake_ 3.4.2
* BWA_ 0.7.12
* SAMtools_ 1.3.1
* BCFtools_ 1.3.1
* Graphviz_ 2.38.0
* PyYAML_ 3.11
* Docutils_ 0.12

The easiest way to setup these prerequisites is to use the Miniconda_ Python 3 distribution.
The tutorial assumes that you are using either Linux or MacOS X.
Both Snakemake and Miniconda work also under Windows, but the Windows shell is too different to be able to provide generic examples.

Setup a Linux VM with Vagrant under Windows
:::::::::::::::::::::::::::::::::::::::::::

If you already use Linux or MacOS X, go on with **Step 1**.
If you use Windows, you can setup a Linux virtual machine (VM) with Vagrant_.
First, install Vagrant following the installation instructions in the `Vagrant Documentation`_.
Then, create a reasonable new directory you want to share with your Linux VM, e.g., create a folder ``vagrant-linux`` somewhere.
Open a command line prompt, and change into that directory.
Here, you create a 64-bit Ubuntu Linux environment with

.. code:: bash

    vagrant init hashicorp/precise64
    vagrant up

If you decide to use a 32-bit image, you will need to download the 32-bit version of Miniconda in the next step.
The contents of the ``vagrant-linux`` folder will be shared with the virtual machine that is set up by vagrant.
You can log into the virtual machine via

.. code:: bash

    vagrant ssh

If this command tells you to install an SSH client, you can follow the instructions in this Blogpost_.
Now, you can follow the steps of our tutorial from within your Linux VM.


Step 1: Installing Miniconda 3
::::::::::::::::::::::::::::::

First, please **open a terminal** or make sure you are logged into your Vagrant Linux VM.
Assuming that you have a 64-bit system, on Linux, download and install Miniconda 3 with

.. code:: bash

    wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
    bash Miniconda3-latest-Linux-x86_64.sh

On MacOS X, download and install with

.. code:: bash

    curl https://repo.continuum.io/miniconda/Miniconda3-latest-MacOSX-x86_64.sh -o Miniconda3-latest-MacOSX-x86_64.sh
    bash Miniconda3-latest-MacOSX-x86_64.sh

For a 32-bit system, URLs and file names are analogous but without the ``_64``.
When you are asked the question

.. code::

    Do you wish the installer to prepend the Miniconda3 install location to PATH ...? [yes|no]

answer with **yes**.
Along with a minimal Python 3 environment, Miniconda contains the package manager Conda_.
After opening a **new terminal**, you can use the new ``conda`` command to install software packages and create isolated environments to, e.g., use different versions of the same package.
We will later use Conda_ to create an isolated enviroment with all required software for this tutorial.

Step 2: Preparing a working directory
:::::::::::::::::::::::::::::::::::::

First, **create a new directory** ``snakemake-tutorial`` at a reasonable place and **change into that directory** in your terminal.
If you use a Vagrant Linux VM from Windows as described above, create the directory under ``/vagrant/``, so that the contents are shared with your host system (you can then edit all files from within Windows with an editor that supports Unix line breaks).
In this directory, we will later create an example workflow that illustrates the Snakemake syntax and execution environment.
First, we download some example data on which the workflow shall be executed:

.. code:: bash

    wget https://bitbucket.org/snakemake/snakemake/downloads/snakemake-tutorial-data.tar.gz
    tar -xf snakemake-tutorial-data.tar.gz

This will create a ``data`` folder and a ``requirements.txt`` file in the working directory.

Step 3: Creating an environment with the required software
::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

The ``requirements.txt`` file can be used to install all required software into an isolated conda environment with the name ``snakemake-tutorial`` via

.. code:: bash

    conda create -n snakemake-tutorial -c bioconda --file requirements.txt

Note that the arguments after the ``-c`` flags define software channels that shall be used in addition to the main ``conda`` repository.
Here, we use the Bioconda_ channel, which contains a growing collection of bioinformatics software packaged for Conda.

Step 4: Activating the environment
::::::::::::::::::::::::::::::::::

To activate the ``snakemake-tutorial`` enviroment, execute

.. code:: bash

    source activate snakemake-tutorial

Now you can use the installed tools.
Execute

.. code:: bash

    snakemake --help

to test this and get information about the command-line interface of Snakemake.
To exit the environment, you can execute

.. code:: bash

    source deactivate

but **don't do that now**, since we finally want to start working with Snakemake :-).

------------

Basics: An example workflow
---------------------------

Please make sure that you have **activated** the environment we created before, and that you have an open terminal in the working directory you have created.

**A Snakemake workflow is defined by specifying rules in a Snakefile**.
**Rules decompose the workflow into small steps** (e.g., the application of a single tool) by specifying how to create sets of **output files** from sets of **input files**.
Snakemake automatically **determines the dependencies** between the rules by matching file names.

The Snakemake language extends the Python language, adding syntactic structures for rule definition and additional controls.
All added syntactic structures begin with a keyword followed by a code block that is either in the same line or indented and consisting of multiple lines.
The resulting syntax resembles that of original Python constructs.

In the following, we will introduce the Snakemake syntax by creating an example workflow.
The workflow will map sequencing reads to a reference genome and call variants on the mapped reads.

Step 1: Mapping reads
:::::::::::::::::::::

Our first Snakemake rule maps reads of a given sample to a given reference genome.
In the working directory, **create a new file** called ``Snakefile`` with an editor of your choice.
We propose to use the Atom_ editor, since it provides out-of-the-box syntax highlighting for Snakemake.
In the Snakefile, define the following rule:

.. code:: bash

    rule bwa_map:
        input:
            "data/genome.fa",
            "data/samples/A.fastq"
        output:
            "mapped_reads/A.bam"
        shell:
            "bwa mem {input} | samtools view -Sb - > {output}"

A Snakemake rule has a name (here ``bwa_map``) and a number of directives, here ``input``, ``output`` and ``shell``.
The ``input`` and ``output`` directives are followed by lists of files that are expected to be used or created by the rule.
In the simplest case, these are just explicit Python strings.
The ``shell`` directive is followed by a Python string containing the shell command to execute.
In the shell command string, we can refer to elements of the rule via braces notation (similar to the Python format function).
Here, we refer to the output file by specifying ``{output}`` and to the input files by specifying ``{input}``.
Since the rule has multiple input files, Snakemake will concatenate them separated by a whitespace.
In other words, Snakemake will replace ``{input}`` with ``data/genome.fa data/samples/A.fastq`` before executing the command.
The shell command invokes ``bwa mem`` with reference genome and reads, and pipes the output into ``samtools`` which creates a compressed BAM file containing the alignments.
The output of ``samtools`` is piped into the output file defined by the rule.

When a workflow is executed, Snakemake tries to generate given **target** files.
Target files can be specified via the command line.
By executing

.. code:: bash

    snakemake -np mapped_reads/A.bam

in the working directory containing the Snakefile, we tell Snakemake to generate the target file ``mapped_reads/A.bam``.
Since we used the ``-n`` (or ``--dryrun``) flag, Snakemake will only show the execution plan instead of actually perform the steps.
The ``-p`` flag instructs Snakemake to also print the resulting shell command for illustation.
To generate the target files, **Snakemake applies the rules given in the Snakefile in a top-down way**.
The application of a rule to generate a set of output files is called **job**.
For each input file of a job, Snakemake again (i.e. recursively) determines rules that can be applied to generate it.
This yields a directed acyclic graph (DAG) of jobs where the edges represent dependencies.
So far, we only have a single rule, and the DAG of jobs consists of a single node.
Nevertheless, we can **execute our workflow** with

.. code:: bash

    snakemake mapped_reads/A.bam

Note that, after completion of above command, Snakemake will not try to create ``mapped_reads/A.bam`` again, because it is already present in the file system.
Snakemake **only re-runs jobs if one of the input files is newer than one of the output files or one of the input files will be updated by another job**.

Step 2: Generalizing the read mapping rule
::::::::::::::::::::::::::::::::::::::::::

Obviously, the rule will only work for a single sample with reads in the file ``data/samples/A.fastq``.
However, Snakemake allows to **generalize rules by using named wildcards**.
Simply replace the ``A`` in the second input file and in the output file with the wildcard ``{sample}``, leading to

.. code:: bash

    rule bwa_map:
        input:
            "data/genome.fa",
            "data/samples/{sample}.fastq"
        output:
            "mapped_reads/{sample}.bam"
        shell:
            "bwa mem {input} | samtools view -Sb - > {output}"

When Snakemake determines that this rule can be applied to generate a target file by replacing the wildcard ``{sample}`` in the output file with an appropriate value, it will propagate that value to all occurences of ``{sample}`` in the input files and thereby determine the necessary input for the resulting job.
Note that you can have multiple wildcards in your file paths, however, to avoid conflicts with other jobs of the same rule, **all output files** of a rule have to **contain exactly the same wildcards**.

When executing

.. code:: bash

    snakemake -np mapped_reads/B.bam

Snakemake will determine that the rule ``bwa_map`` can be applied to generate the target file by replacing the wildcard ``{sample}`` with the value ``B``.
In the output of the dry-run, you will see how the wildcard value is propagated to the input files and all filenames in the shell command.
You can also **specify multiple targets**, e.g.:

.. code:: bash

    snakemake -np mapped_reads/A.bam mapped_reads/B.bam

Some Bash_ magic can make this particularly handy. For example, you can alternatively compose our multiple targets in a single pass via

.. code:: bash

    snakemake -np mapped_reads/{A,B}.bam

Note that this is not a special Snakemake syntax. Bash is just expanding the given path into two, one for each element of the set ``{A,B}``.

In both cases, you will see that Snakemake only proposes to create the output file ``mapped_reads/B.bam``.
This is because you already executed the workflow before (see the previous step) and no input file is newer than the output file ``mapped_reads/A.bam``.
You can update the file modification date of the input file
``data/samples/A.fastq`` via

.. code:: bash

    touch data/samples/A.fastq

and see how Snakemake wants to re-run the job to create the file ``mapped_reads/A.bam`` by executing

.. code:: bash

    snakemake -np mapped_reads/A.bam mapped_reads/B.bam


Step 3: Sorting read alignments
:::::::::::::::::::::::::::::::

For later steps, we need the read alignments in the BAM files to be sorted.
This can be achieved with the ``samtools`` command.
We add the following rule beneath the ``bwa_map`` rule:

.. code:: bash

    rule samtools_sort:
        input:
            "mapped_reads/{sample}.bam"
        output:
            "sorted_reads/{sample}.bam"
        shell:
            "samtools sort -T sorted_reads/{wildcards.sample} "
            "-O bam {input} > {output}"

This rule will take the input file from the ``mapped_reads`` directory and store a sorted version in the ``sorted_reads`` directory.
Note that Snakemake **automatically creates missing directories** before jobs are executed.
For sorting, ``samtools`` requires a prefix specified with the flag ``-T``.
Here, we need the value of the wildcard ``sample``.
Snakemake allows to access wildcards in the shell command via the ``wildcards`` object that has an attribute with the value for each wildcard.

When issuing

.. code:: bash

    snakemake -np sorted_reads/B.bam

you will see how Snakemake wants to run first the rule ``bwa_map`` and then the rule ``samtools_sort`` to create the desired target file:
as mentioned before, the dependencies are resolved automatically by matching file names.

Step 4: Indexing read alignments and visualizing the DAG of jobs
::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

Next, we need to index the sorted read alignments for random access.
This can be done with the following rule:

.. code:: bash

    rule samtools_index:
        input:
            "sorted_reads/{sample}.bam"
        output:
            "sorted_reads/{sample}.bam.bai"
        shell:
            "samtools index {input}"

Having three steps already, it is a good time to take a closer look at the resulting DAG of jobs.
By executing

.. code:: bash

    snakemake --dag sorted_reads/{A,B}.bam.bai | dot -Tsvg > dag.svg

we create a **visualization of the DAG** using the ``dot`` command provided by Graphviz_.
For the given target files, Snakemake specifies the DAG in the dot language and pipes it into the ``dot`` command, which renders the definition into SVG format.
The rendered DAG is piped into the file ``dag.svg`` and will look similar to this:

.. image::workflow/dag_index.png
   :align: center

The DAG contains a node for each job and edges representing the dependencies.
Jobs that don't need to be run because their output is up-to-date are dashed.
For rules with wildcards, the value of the wildcard for the particular job is displayed in the job node.

Exercise
........

* Run parts of the workflow using different targets. Recreate the DAG and see how different rules become dashed because their output is present and up-to-date.

Step 5: Calling genomic variants
::::::::::::::::::::::::::::::::

The next step in our workflow will aggregate the aligned reads from all samples and jointly call genomic variants on them.
Snakemake provides a **helper function for collecting input files**.
With

.. code:: bash

    expand("sorted_reads/{sample}.bam", sample=SAMPLES)

we obtain a list of files where the given pattern ``"sorted_reads/{sample}.bam"`` was formatted with the values in the given list of samples ``SAMPLES``, i.e.

.. code:: bash

    ["sorted_reads/A.bam", "sorted_reads/B.bam"]

The function is particularly useful when the pattern contains multiple wildcards.
For example,

.. code:: bash

    expand("sorted_reads/{sample}.{replicate}.bam", sample=SAMPLES, replicate=[0, 1])

would create the product of all elements of ``SAMPLES`` and the list ``[0, 1]``, yielding

.. code:: bash

    ["sorted_reads/A.0.bam", "sorted_reads/A.1.bam", "sorted_reads/B.0.bam", "sorted_reads/B.1.bam"]

For more information, see the Documentation_.
Here, we use only the simple case of ``expand``.
We first let Snakemake know which samples we want to consider.
Remember that Snakemake works top-down, it does not automatically infer this from, e.g., the fastq files in the data folder.
Remember that Snakefiles are in principle Python code enhanced by some declarative statements to define workflows.
Hence, we can define the list of samples ad-hoc in plain Python at the top of the Snakefile:

.. code:: bash

    SAMPLES = ["A", "B"]

Later, we will learn about more sophisticated ways like **config files**.
Now, we can add the following rule to our Snakefile:

.. code:: bash

    rule bcftools_call:
        input:
            fa="data/genome.fa",
            bam=expand("sorted_reads/{sample}.bam", sample=SAMPLES),
            bai=expand("sorted_reads/{sample}.bam.bai", sample=SAMPLES)
        output:
            "calls/all.vcf"
        shell:
            "samtools mpileup -g -f {input.fa} {input.bam} | "
            "bcftools call -mv - > {output}"

With multiple input or output files, it is sometimes handy to refer them separately in the shell command.
This can be done by **specifying names for input or output files** (here, e.g., ``fa=...``).
The files can then be referred in the shell command via, e.g., ``{input.fa}``.
For **long shell commands** like this one, it is advisable to **split the string over multiple indented lines**.
Python will automatically merge it into one.
Further, you will notice that the **input or output file lists can contain arbitrary Python statements**, as long as it returns a string, or a list of strings.
Here, we invoke our ``expand`` function to aggregate over the aligned reads of all samples.

Exercise
........

* obtain the updated DAG of jobs for the target file ``calls/all.vcf``, it should look like this:

.. image::workflow/dag_call.png
   :align: center

Step 6: Writing a report
::::::::::::::::::::::::

Although Snakemake workflows are already self-documenting to a certain degree, it is often useful to summarize the obtained results and performed steps in a comprehensive **report**.
With Snakemake, such reports can be composed easily with the built-in ``report`` function.
It is best practice to create reports in a separate rule that takes all desired results as input files and provides a **single HTML file as output**.

.. code:: bash

    rule report:
        input:
            "calls/all.vcf"
        output:
            "report.html"
        run:
            from snakemake.utils import report
            with open(input[0]) as vcf:
                n_calls = sum(1 for l in vcf if not l.startswith("#"))

            report("""
            An example variant calling workflow
            ===================================

            Reads were mapped to the Yeast
            reference genome and variants were called jointly with
            SAMtools/BCFtools.

            This resulted in {n_calls} variants (see Table T1_).
            """, output[0], T1=input[0])

First, we notice that this rule does not entail a shell command.
Instead, we use the ``run`` directive, which is followed by plain Python code.
Similar to the shell case, we have access to ``input`` and ``output`` files, which we can handle as plain Python objects (no braces notation here).

We go through the ``run`` block line by line.
First, we import the ``report`` function from ``snakemake.utils``.
Second, we open the VCF file by accessing it via its index in the input files (i.e. ``input[0]``), and count the number of non-header lines (which is equivalent to the number of variant calls).
Third, we create the report using the ``report`` function.
The function takes a string that contains RestructuredText_ markup.
In addition, we can use the familiar braces notation to access any Python variables (here the ``samples`` and ``n_calls`` variables we have defined before).
The second argument of the ``report`` function is the path were the report will be stored (the function creates a single HTML file).
Then, report expects any number of keyword arguments referring to files that shall be embedded into the report.
Technically, this means that the file will be stored as a Base64 encoded `data URI`_ within the HTML file, making reports entirely self-contained.
Importantly, you can refer to the files from within the report via the given keywords followed by an underscore (here ``T1_``).
Hence, reports can be used to semantically connect and explain the obtained results.

When having many result files, it is sometimes handy to define the names already in the list of input files and unpack these into keyword arguments as follows:

.. code:: bash

    report("""...""", output[0], **input)

Further, you can add meta data in the form of any string that will be displayed in the footer of the report, e.g.

.. code:: bash

    report("""...""", output[0], metadata="Author: Johannes Köster (koester@jimmy.harvard.edu)", **input)


Step 7: Adding a target rule
::::::::::::::::::::::::::::

So far, we always executed the workflow by specifying a target file at the command line.
Apart from filenames, Snakemake **also accepts rule names as targets** if the referred rule does not have wildcards.
Hence, it is possible to write target rules collecting particular subsets of the desired results or all results.
Moreover, if no target is given at the command line, Snakemake will define the **first rule** of the Snakefile as the target.
Hence, it is best practice to have a rule ``all`` at the top of the workflow which has all typically desired target files as input files.

Here, this means that we add a rule

.. code:: bash

    rule all:
        input:
            "report.html"

to the top of our workflow.
When executing Snakemake with

.. code:: bash

    snakemake -n

the execution plan for creating the file ``report.html`` which contains and summarizes all our results will be shown.
Note that, apart from Snakemake considering the first rule of the workflow as default target, **the appearance of rules in the Snakefile is arbitrary and does not influence the DAG of jobs**.

Exercise
........

* Create the DAG of jobs for the complete workflow.
* Execute the complete workflow and have a look at the resulting ``report.html`` in your browser.
* Snakemake provides handy flags for forcing re-execution of parts of the workflow. Have a look at the command line help with ``snakemake --help`` and search for the flag ``--forcerun``. Then, use this flag to re-execute the rule ``samtools_sort`` and see what happens.
* With ``--reason`` it is possible to display the execution reason for each job. Try this flag together with a dry-run and the ``--forcerun`` flag to understand the decisions of Snakemake.

Summary
:::::::

In total, the resulting workflow looks like this:

.. code:: bash

    SAMPLES = ["A", "B"]


    rule all:
        input:
            "report.html"


    rule bwa_map:
        input:
            "data/genome.fa",
            "data/samples/{sample}.fastq"
        output:
            "mapped_reads/{sample}.bam"
        shell:
            "bwa mem {input} | samtools view -Sb - > {output}"


    rule samtools_sort:
        input:
            "mapped_reads/{sample}.bam"
        output:
            "sorted_reads/{sample}.bam"
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
            bam=expand("sorted_reads/{sample}.bam", sample=SAMPLES),
            bai=expand("sorted_reads/{sample}.bam.bai", sample=SAMPLES)
        output:
            "calls/all.vcf"
        shell:
            "samtools mpileup -g -f {input.fa} {input.bam} | "
            "bcftools call -mv - > {output}"


    rule report:
        input:
            "calls/all.vcf"
        output:
            "report.html"
        run:
            from snakemake.utils import report
            with open(input[0]) as vcf:
                n_calls = sum(1 for l in vcf if not l.startswith("#"))

            report("""
            An example variant calling workflow
            ===================================

            Reads were mapped to the Yeast
            reference genome and variants were called jointly with
            SAMtools/BCFtools.

            This resulted in {n_calls} variants (see Table T1_).
            """, output[0], T1=input[0])


------------

Advanced: Decorating the example workflow
-----------------------------------------

Now that the basic concepts of Snakemake have been illustrated, we can introduce advanced topics.

Step 1: Specifying the number of used threads
:::::::::::::::::::::::::::::::::::::::::::::

For some tools, it is advisable to use more than one thread in order to speed up the computation.
**Snakemake can be made aware of the threads a rule needs** with the ``threads`` directive.
In our example workflow, it makes sense to use multiple threads for the rule ``bwa_map``:

.. code:: bash

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

.. code:: bash

    snakemake --cores 10

would execute the workflow with 10 cores.
Since the rule ``bwa_map`` needs 8 threads, only one job of the rule can run at a time, and the Snakemake scheduler will try to saturate the remaining cores with other jobs like, e.g., ``samtools_sort``.
The threads directive in a rule is interpreted as a maximum: when **less cores than threads** are provided, the number of threads a rule uses will be **reduced to the number of given cores**.

Apart from the very common thread resource, Snakemake provides a ``resources`` directive that can be used to **specify arbitrary resources**, e.g., memory usage or auxiliary computing devices like GPUs.
Similar to threads, these can be considered by the scheduler when an available amount of that resource is given with the command line argument ``--resources``.
Details can be found in the Snakemake Documentation_.

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

.. code:: bash

    configfile: "config.yaml"

to the top of the Snakefile.
Snakemake will load the config file and store its contents into a globally available dictionary named ``config``.
In our case, it makes sense to specify the samples in ``config.yaml`` as

.. code:: yaml

    samples:
        A: data/samples/A.fastq
        B: data/samples/B.fastq

Now, we can remove the statement defining ``SAMPLES`` from the Snakefile and change the rule ``bcftools_call`` to

.. code:: bash

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

.. code:: bash

    rule bwa_map:
        input:
            "data/genome.fa",
            lambda wildcards: config["samples"][wildcards.sample]
        output:
            "mapped_reads/{sample}.bam"
        threads: 8
        shell:
            "bwa mem -t {threads} {input} | samtools view -Sb - > {output}"

Here, we use an anonymous function, also called **lambda expression**.
Any normal function would work as well.
Input functions take as **single argument** a ``wildcards`` object, that allows to access the wildcards values via attributes (here ``wildcards.sample``).
They **return a string or a list of strings**, that are interpeted as paths to input files (here, we return the path that is stored for the sample in the config file).
Input functions are evaluated once the wildcard values of a job are determined.


Exercise
........

* In the ``data/samples`` folder, there is an additional sample ``C.fastq``. Add that sample to the config file and see how Snakemake wants to recompute the part of the workflow belonging to the new sample.


Step 4: Rule parameters
:::::::::::::::::::::::

Sometimes, shell commands are not only composed of input and output files and some static flags.
In particular, it can happen that additional parameters need to be set depending on the wildcard values of the job.
For this, Snakemake allows to **define arbitrary parameters** for rules with the ``params`` directive.
In our workflow, it is reasonable to annotate aligned reads with so-called read groups, that contain metadata like the sample name.
We modify the rule ``bwa_map`` accordingly:

.. code:: bash

    rule bwa_map:
        input:
            "data/genome.fa",
            lambda wildcards: config["samples"][wildcards.sample]
        output:
            "mapped_reads/{sample}.bam"
        params:
            rg="@RG\tID:{sample}\tSM:{sample}"
        threads: 8
        shell:
            "bwa mem -R '{params.rg}' -t {threads} {input} | samtools view -Sb - > {output}"

Similar to input and output files, ``params`` can be accessed from the shell command.
Moreover, the ``params`` directive can also take functions like in Step 3 to defer initialization to the DAG phase.

Exercise
........

* Variant calling can consider a lot of parameters. A particularly important one is the prior mutation rate (1e-3 per default). It is set via the flag ``-P`` of the ``bcftools call`` command. Consider making this flag configurable via adding a new key to the config file and using the ``params`` directive in the rule ``bcftools_call`` to propagate it to the shell command.

Step 5: Logging
:::::::::::::::

When executing a large workflow, it is usually desirable to store the output of each job persistently in files instead of just printing it to the terminal.
For this purpose, Snakemake allows to **specify log files** for rules.
Log files are defined via the ``log`` directive and handled similarly to output files, but they are not subject of rule matching and are not cleaned up when a job fails.
We modify our rule ``bwa_map`` as follows:

.. code:: bash

    rule bwa_map:
        input:
            "data/genome.fa",
            lambda wildcards: config["samples"][wildcards.sample]
        output:
            "mapped_reads/{sample}.bam"
        params:
            rg="@RG\tID:{sample}\tSM:{sample}"
        log:
            "logs/bwa_map/{sample}.log"
        threads: 8
        shell:
            "(bwa mem -R '{params.rg}' -t {threads} {input} | "
            "samtools view -Sb - > {output}) 2> {log}"

The shell command is modified to collect STDERR output of both ``bwa`` and ``samtools`` and pipe it into the file referred by ``{log}``.
It is best practice to store all log files in a ``logs`` subdirectory, prefixed by the rule or tool name.
Log files must contain exactly the same wildcards as the output files to avoid clashes.

Exercise
........

* Add a log directive to the ``bcftools_call`` rule as well.
* Time to re-run the whole workflow (remember the command line flags to force re-execution). See how log files are created for variant calling and read mapping.
* The ability to track the provenance of each generated result is an important step towards reproducible analyses. Apart from the ``report`` functionality discussed before, Snakemake can summarize various provenance information for all output files of the workflow. The flag ``--summary`` prints a table associating each output file with the rule used to generate it, the creation date and optionally the version of the tool used for creation is provided. Further, the table informs about updated input files and changes to the source code of the rule after creation of the output file. Invoke Snakemake with ``--summary`` to examine the information for our example.

Step 6: Temporary and protected files
:::::::::::::::::::::::::::::::::::::

In our workflow, we create two BAM files for each sample, namely
the output of the rules ``bwa_map`` and ``samtools_sort``.
When not dealing with examples, the underlying data is usually huge.
Hence, the resulting BAM files need a lot of disk space and their creation takes some time.
Snakemake allows to **mark output files as temporary**, such that they are deleted once every consuming job has been executed, in order to save disk space.
We use this mechanism for the output file of the rule ``bwa_map``:

.. code:: bash

    rule bwa_map:
        input:
            "data/genome.fa",
            lambda wildcards: config["samples"][wildcards.sample]
        output:
            temp("mapped_reads/{sample}.bam")
        params:
            rg="@RG\tID:{sample}\tSM:{sample}"
        log:
            "logs/bwa_map/{sample}.log"
        threads: 8
        shell:
            "(bwa mem -R '{params.rg}' -t {threads} {input} | "
            "samtools view -Sb - > {output}) 2> {log}"

This results in the deletion of the BAM file once the corresponding ``samtools_sort`` job has been executed.
Since the creation of BAM files via read mapping and sorting is computationally expensive, it is reasonable to **protect** the final BAM file **from accidental deletion or modification**.
We modify the rule ``samtools_sort`` by marking it's output file as ``protected``:

.. code:: bash

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

.. code:: bash

    configfile: "config.yaml"


    rule all:
        input:
            "report.html"


    rule bwa_map:
        input:
            "data/genome.fa",
            lambda wildcards: config["samples"][wildcards.sample]
        output:
            temp("mapped_reads/{sample}.bam")
        params:
            rg="@RG\tID:{sample}\tSM:{sample}"
        log:
            "logs/bwa_map/{sample}.log"
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


    rule report:
        input:
            "calls/all.vcf"
        output:
            "report.html"
        run:
            from snakemake.utils import report
            with open(input[0]) as vcf:
                n_calls = sum(1 for l in vcf if not l.startswith("#"))

            report("""
            An example variant calling workflow
            ===================================

            Reads were mapped to the Yeast
            reference genome and variants were called jointly with
            SAMtools/BCFtools.

            This resulted in {n_calls} variants (see Table T1_).
            """, output[0], T1=input[0])

------------

Additional features
-------------------

In the following, we introduce some features that are beyond the scope of above example workflow.
For details and even more features, see the Documentation_, the FAQ_ and the command line help (``snakemake --help``).


Benchmarking
::::::::::::

With the ``benchmark`` directive, Snakemake can be instructed to **measure the wall clock time of a job**.
We activate benchmarking for the rule ``bwa_map``:

.. code:: bash

    rule bwa_map:
        input:
            "data/genome.fa",
            lambda wildcards: config["samples"][wildcards.sample]
        output:
            temp("mapped_reads/{sample}.bam")
        params:
            rg="@RG\tID:{sample}\tSM:{sample}"
        log:
            "logs/bwa_map/{sample}.log"
        benchmark:
            "benchmarks/{sample}.bwa.benchmark.txt"
        threads: 8
        shell:
            "(bwa mem -R '{params.rg}' -t {threads} {input} | "
            "samtools view -Sb - > {output}) 2> {log}"

The ``benchmark`` directive takes a string that points to the file where benchmarking results shall be stored.
Similar to output files, the path can contain wildcards (it must be the same wildcards as in the output files).
When a job derived from the rule is executed, Snakemake will measure the wall clock time and store it in the file in tab-delimited format.
With the command line flag ``--benchmark-repeats``, Snakemake can be instructed to perform repetetive measurements by executing benchmark jobs multiple times.
The repeated measurements occur as subsequent lines in the tab-delimited benchmark file.

We can include the benchmark results into our report:

.. code:: bash

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

.. code:: bash

    include: "path/to/other.snakefile"

Alternatively, Snakemake allows to **define sub-workflows**.
A sub-workflow refers to a working directory with a complete Snakemake workflow.
Output files of that sub-workflow can be used in the current Snakefile.
When executing, Snakemake ensures that the output files of the sub-workflow are up-to-date before executing the current workflow.
This mechanism is particularly useful when you want to extend a previous analysis without modifying it.
For details about sub-workflows, see the Documentation_.


Exercise
........

* Put the read mapping related rules into a separate Snakefile and use the ``include`` directive to make them available in our example workflow again.


Using custom scripts
::::::::::::::::::::

Using the ``run`` directive as above is only reasonable for short Python scripts.
As soon as your script becomes larger, it is reasonable to separate it from the
workflow definition.
For this purpose, Snakemake offers the ``script`` directive.
Using this, ``report`` rule from above could instead look like this:

.. code:: bash

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

For details and examples, see the `script section in the Documentation`_.

Cluster execution
:::::::::::::::::

By default, Snakemake executes jobs on the local machine it is invoked on.
Alternatively, it can execute jobs in **distributed environments, e.g., compute clusters or batch systems**.
If the nodes share a common file system, Snakemake supports three alternative execution modes.

In cluster enviroments, compute jobs are usually submitted as shell scripts via commands like ``qsub``.
Snakemake provides a **generic mode** to execute on such clusters.
By invoking Snakemake with

.. code:: bash

    snakemake --cluster qsub --jobs 100

each job will be compiled into a shell script that is submitted with the given command (here ``qsub``).
The ``--jobs`` flag limits the number of concurrently submitted jobs to 100.
This basic mode assumes that the submission command returns immediately after submitting the job.
Some clusters allow to run the submission command in **synchronous mode**, such that it waits until the job has been executed.
In such cases, we can invoke e.g.

.. code:: bash

    snakemake --cluster-sync "qsub -sync yes" --jobs 100

The specified submission command can also be **decorated with additional parameters taken from the submitted job**.
For example, the number of used threads can be accessed in braces similarly to the formatting of shell commands, e.g.

.. code:: bash

    snakemake --cluster "qsub -pe threaded {threads}" --jobs 100

Alternatively, Snakemake can use the Distributed Resource Management Application API (DRMAA_).
This API provides a common interface to control various resource management systems.
The **DRMAA support** can be activated by invoking Snakemake as follows:

.. code:: bash

    snakemake --drmaa --jobs 100

If available, **DRMAA is preferable over the generic cluster modes** because it provides better control and error handling.
To support additional cluster specific parametrization, a Snakefile can be complementd by a `cluster configuration`_.


Constraining wildcards
::::::::::::::::::::::

Snakemake uses regular expressions to match output files to input files and determine dependencies between the jobs.
Sometimes it is useful to constrain the values a wildcard can have.
This can be achieved by adding a regular expression that describes the set of allowed wildcard values.
For example, the wildcard ``sample`` in the output file ``"sorted_reads/{sample}.bam"`` can be constrained to only allow alphanumeric sample names as ``"sorted_reads/{sample,[A-Za-z0-9]+}.bam"``.
This mechanism helps to solve two kinds of ambiguity.

* It can help to avoid ambiguous rules, i.e. two or more rules that can be applied to generate the same output file. Other ways of handling ambiguous rules are described in the Documentation_.
* It can help to guide the regular expression based matching so that wildcards are assigned to the right parts of a file name. Consider the output file ``{sample}.{group}.txt`` and assume that the target file is ``A.1.normal.txt``. It is not clear whether ``dataset="A.1"`` and ``group="normal"`` or ``dataset="A"`` and ``group="1.normal"`` is the right assignment. Here, constraining the dataset wildcard by ``{sample,[A-Z]+}.{group}`` solves the problem.

When dealing with ambiguous rules, it is best practice to first try to solve the ambiguity by using a proper file structure, for example, by separating the output files of different steps in different directories.
