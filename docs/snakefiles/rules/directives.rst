
.. _snakefiles-directives:

Directives
----------

.. _snakefiles-threads:

``threads``
^^^^^^^^^^

Further, a rule can be given a number of threads to use, i.e.

.. code-block:: python

    rule NAME:
        input: "path/to/inputfile", "path/to/other/inputfile"
        output: "path/to/outputfile", "path/to/another/outputfile"
        threads: 8
        shell: "somecommand --threads {threads} {input} {output}"

.. note::

    On a cluster node, Snakemake uses as many cores as available on that node.
    Hence, the number of threads used by a rule never exceeds the number of physically available cores on the node.
    Note: This behavior is not affected by ``--local-cores``, which only applies to jobs running on the main node.

Snakemake can alter the number of cores available based on command line options. Therefore it is useful to propagate it via the built in variable ``threads`` rather than hardcoding it into the shell command.
In particular, it should be noted that the specified threads have to be seen as a maximum. When Snakemake is executed with fewer cores, the number of threads will be adjusted, i.e. ``threads = min(threads, cores)`` with ``cores`` being the number of cores specified at the command line (option ``--cores``).

Hardcoding a particular maximum number of threads like above is useful when a certain tool has a natural maximum beyond which parallelization won't help to further speed it up.
This is often the case, and should be evaluated carefully for production workflows.
Also, setting a ``threads:`` maximum is required to achieve parallelism in tools that (often implicitly and without the user knowing) rely on an environment variable for the maximum of cores to use.
For example, this is the case for many linear algebra libraries and for OpenMP.
Snakemake limits the respective environment variables to one core by default, to avoid unexpected and unlimited core-grabbing, but will override this with the ``threads:`` you specify in a rule (the parameters set to ``threads:``, or defaulting to ``1``, are: ``OMP_NUM_THREADS``, ``GOTO_NUM_THREADS``, ``OPENBLAS_NUM_THREADS``, ``MKL_NUM_THREADS``, ``VECLIB_MAXIMUM_THREADS``, ``NUMEXPR_NUM_THREADS``).

If it is certain that no maximum for efficient parallelism exists for a tool, one can instead define threads as a function of the number of cores given to Snakemake:

.. code-block:: python

    rule NAME:
        input: "path/to/inputfile", "path/to/other/inputfile"
        output: "path/to/outputfile", "path/to/another/outputfile"
        threads: workflow.cores * 0.75
        shell: "somecommand --threads {threads} {input} {output}"

The number of given cores is globally available in the Snakefile as an attribute of the workflow object: ``workflow.cores``.
Any arithmetic operation can be performed to derive a number of threads from this. E.g., in the above example, we reserve 75% of the given cores for the rule.
Snakemake will always round the calculated value down (while enforcing a minimum of 1 thread).

Starting from version 3.7, threads can also be a callable that returns an ``int`` value. The signature of the callable should be ``callable(wildcards[, input])`` (input is an optional parameter).  It is also possible to refer to a predefined variable (e.g, ``threads: threads_max``) so that the number of cores for a set of rules can be changed with one change only by altering the value of the variable ``threads_max``.

Both threads can be defined (or overwritten) upon invocation (without modifying the workflow code) via `--set-threads` see :ref:`all_options` and via workflow profiles, see :ref:`executing-profiles`.
To quickly exemplify the latter, you could provide the following workflow profile in a file ``profiles/default/config.yaml`` relative to the Snakefile or the current working directory:

.. code-block:: yaml

    set-threads:
        b: 16

to set the (maximum) number of threads rule ``b`` uses to 16.



.. _snakefiles-resources:

``resources``
^^^^^^^^^^^^^

In addition to threads, a rule can use arbitrary user-defined resources by specifying them with the resources-keyword:

.. code-block:: python

    rule a:
        input:     ...
        output:    ...
        resources:
            mem_mb=100
        shell:
            "..."

If workflow-wide limits for the resources are given via the command line, e.g.

.. code-block:: console

    $ snakemake --resources mem_mb=200


the scheduler will ensure that the given resources are not exceeded by running jobs.
Resources are always meant to be specified as total per job, not by thread (i.e. above ``mem_mb=100`` in rule ``a`` means that any job from rule ``a`` will require ``100`` megabytes of memory in total, and not per thread).

**Importantly**, there are some :ref:`standard resources <snakefiles-standard-resources>` that should be considered before making up your own.

In general, resources are just names to the Snakemake scheduler, i.e., Snakemake does not check on the resource consumption of jobs in real time.
Instead, resources are used to determine which jobs can be executed at the same time without exceeding the limits specified at the command line.
Apart from making Snakemake aware of hybrid-computing architectures (e.g. with a limited number of additional devices like GPUs) this allows us to control scheduling in various ways, e.g. to limit IO-heavy jobs by assigning an artificial IO-resource to them and limiting it via the ``--resources`` flag.
If no limits are given, the resources are ignored in local execution.

Resources can have any arbitrary name, and must be assigned ``int`` or ``str`` values.
In case of ``None``, the resource is considered to be unset (i.e. ignored) in the rule.

.. _snakefiles-dynamic-resources:

Dynamic Resources
"""""""""""""""""

It is often useful to determine resource specifications dynamically during workflow execution.
A common example is determining the amount of memory that a job needs, based on the input file size of that particular rule instance.
To enable this, resource specifications can also be callables (for example functions or lambda expressions) that return ``int``, ``str`` or ``None`` values.
The signature of the callable must be ``callable(wildcards [, input] [, threads] [, attempt])`` (``input``, ``threads``, and ``attempt`` are optional parameters).
Such callables are evaluated immediately before the job is executed (or printed during a dry-run).

The above described example of using input size to determined memory requirements could for example be realized via a lambda expression (here also providing a minimum value of 300 MB memory):

.. code-block:: python

    rule:
        input:    ...
        output:   ...
        resources:
            mem_mb=lambda wc, input: max(2.5 * input.size_mb, 300)
        shell:
            "..."

In order to make this work with a dry-run, where the input files are not yet present, Snakemake automatically converts a ``FileNotFoundError`` that is raised by the callable into a placeholder called ``<TBD>`` that will be displayed during dry-run in such a case.

The parameter ``attempt`` allows us to adjust resources based on how often the job has been restarted (see :ref:`all_options`, option ``--retries``).
This is handy when executing a Snakemake workflow in a cluster environment, where jobs can e.g. fail because of too limited resources.
When Snakemake is executed with ``--retries 3``, it will try to restart a failed job 3 times before it gives up.
Thereby, the parameter ``attempt`` will contain the current attempt number (starting from ``1``).
This can be used to adjust the required memory as follows

.. code-block:: python

    def get_mem_mb(wildcards, attempt):
        return attempt * 100

    rule:
        input:    ...
        output:   ...
        resources:
            mem_mb=get_mem_mb
        shell:
            "..."

Here, the first attempt will require 100 MB memory, the second attempt will require 200 MB memory and so on.
When passing memory requirements to the cluster engine, you can by this automatically try out larger nodes if it turns out to be necessary.

Another application of callables as resources is when memory usage depends on the number of threads:

.. code-block:: python

    def get_mem_mb(wildcards, threads):
        return threads * 150

    rule b:
        input:     ...
        output:    ...
        threads: 8
        resources:
            mem_mb=get_mem_mb
        shell:
            "..."

Here, the value that the function ``get_mem_mb`` returns, grows linearly with the number of threads.
Of course, any other arithmetic could be performed in that function.

Both threads and resources can be defined (or overwritten) upon invocation (without modifying the workflow code) via `--set-threads` and `--set-resources`, see :ref:`all_options`.
Or they can be defined via workflow :ref:`executing-profiles`, with the variables listed above in the signature for usable callables.
You could, for example, provide the following workflow profile in a file ``profiles/default/config.yaml`` relative to the Snakefile or the current working directory:

.. code-block:: yaml

    set-threads:
        b: 3
    set-resources:
        b:
            mem_mb: 1000

to set the requirements for rule ``b`` to 3 threads and 1000 MB.

.. _snakefiles-standard-resources:

Standard Resources
""""""""""""""""""

There are several **standard resources**, for total memory, disk usage, runtime, and the temporary directory of a job: ``mem``, ``disk``, ``runtime``, and ``tmpdir``.
All of these resources have specific meanings understood by snakemake and are treated in varying unique ways:

* The ``tmpdir`` resource automatically leads to setting the ``$TMPDIR`` variable for shell commands, scripts, wrappers and notebooks. In cluster or cloud setups, its evaluation is delayed until the actual execution of the job. This way, it can dynamically react on the context of the node of execution.

* The ``runtime`` resource indicates the amount of wall clock time a job needs to run.
  It can be given as string defining a time span or as integer defining **minutes**.
  In the former case, the time span can be defined as a string with a number followed by a unit
  (``ms``, ``s``, ``m``, ``h``, ``d``, ``w``, ``y`` for seconds, minutes, hours, days, and years, respectively).
  The interpretation happens via the `humanfriendly package <https://humanfriendly.readthedocs.io/en/latest/api.html?highlight=parse_timespan#humanfriendly.parse_timespan>`__.
  Cluster or cloud backends may use this to constrain the allowed execution time of the submitted job.
  See :ref:`the section below <resources-remote-execution>` for more information.

* ``disk`` and ``mem`` define the amount of memory and disk space needed by the job.
  They are given as strings with a number followed by a unit (``B``, ``KB``, ``MB``, ``GB``, ``TB``, ``PB``, ``KiB``, ``MiB``, ``GiB``, ``TiB``, ``PiB``).
  The interpretation of the definition happens via the `humanfriendly package <https://humanfriendly.readthedocs.io/en/latest/api.html?highlight=parse_timespan#humanfriendly.parse_size>`__.
  Alternatively, the two can be directly defined as integers via the resources ``mem_mb`` and ``disk_mb`` (to which ``disk`` and ``mem`` are also automatically translated internally).
  They are both locally scoped by default, a fact important for cluster and compute execution.
  :ref:`See below <resources-remote-execution>` for more info.
  They are usually passed to execution backends, e.g. to allow the selection of appropriate compute nodes for the job execution.

* ``gpu``, ``gpu_manufacturer``, and ``gpu_model`` define the number of GPUs, the manufacturer of the GPUs, and the gpu model needed by the job.
  The ``gpu`` resource is an integer and the other two are strings. Please check the executor plugin docs in order to see
  whether and how these resources are supported and properly interpreted by the executor.
  For example, the `kubernetes executor plugin <https://snakemake.github.io/snakemake-plugin-catalog/plugins/executor/kubernetes.html>`__ accepts the terms ``nvidia`` or ``amd`` for the ``gpu_manufacturer`` resource.

Because of these special meanings, the above names should always be used instead of possible synonyms (e.g. ``tmp``, ``time``, ``temp``, etc).

.. _default-resources:

Default Resources
"""""""""""""""""

Since it could be cumbersome to define these standard resources for every rule, you can set default values via the command line flag ``--default-resources`` or in a :ref:`profile <profiles>`.
As with ``--set-resources``, this can be done dynamically, using the variables specified for the callables in the section on :ref:`snakefiles-dynamic-resources`.
If those resource definitions are mandatory for a certain execution mode, Snakemake will fail with a hint if they are missing.
Any resource definitions inside a rule override what has been defined with ``--default-resources``.
If ``--default-resources`` are specified without any further arguments, Snakemake uses ``'mem_mb=min(max(2*input.size_mb, 1000), 8000)'``, ``'disk_mb=max(2*input.size_mb, 1000) if input else 50000'``, and ``'tmpdir=system_tmpdir'``.

* The ``tmpdir`` value points to whatever is the default of the operating system or specified by any of the environment variables ``$TMPDIR``, ``$TEMP``, or ``$TMP`` as outlined `here <https://docs.python.org/3/library/tempfile.html#tempfile.gettempdir>`__.
* The rationale for the default value of ``disk_mb`` is the following: if there are input files, we assume the rule will use at most twice their size during execution.
  If there are no input files, we cannot know what the rule will need, hence we assume a conservative default of 50GB.
* The rationale for the default value of ``mem_mb`` is the following: we try to scale the required memory with the input file size (conservatively assuming that they are loaded entirely into memory).
  However, we stop at 8GB, in order to avoid artificially high requests. Tools that read very large files rather tend to stream them instead of fully loading them into memory.
* If ``--default-resources`` is specified with some definitions, but any of the above defaults (e.g. ``mem_mb``) is omitted, these are still used.
  In order to explicitly unset these defaults, assign them a value of ``None``, e.g. ``--default-resources mem_mb=None``.
* Of course, any rule specifying concrete resources either via the rule definition or via ``--set-resources`` will override the defaults.

.. _resources-remote-execution:

Resources and Remote Execution
""""""""""""""""""""""""""""""

New to Snakemake 7.11. In cluster or cloud execution, resources may represent either a global constraint across all submissions (e.g. number of API calls per second), or a constraint local to each specific job sumbmission (e.g. the amount of memory available on a node).
Snakemake distinguishes between these two types of constraints using **resource scopes**.
By default, ``mem_mb``, ``disk_mb``, and ``threads`` are all considered ``"local"`` resources, meaning specific to individual submissions.
So if a constraint of 16G of memory is given to snakemake (e.g. ``snakemake --resources mem_mb=16000``), each group job will be allowed 16G of memory.
All other resources are considered ``"global"``, meaning they are tracked across all jobs across all submissions.
For example, if ``api_calls`` was limited to 5 and each job scheduled used 1 api call, only 5 jobs would be scheduled at a time, even if more job submissions were available.

These resource scopes may be modified both in the Snakefile and via the CLI parameter ``--set-resource-scopes``.
The CLI parameter takes priority.
Modification in the Snakefile uses the following syntax:

.. code-block:: python

    resource_scopes:
        gpu="local",
        foo="local",
        disk_mb="global"

Here, we set both ``gpu`` and ``foo`` as local resources, and we changed ``disk_mb`` from its default to be a ``global`` resource.
These options could be overridden at the command line using:

.. code-block:: console

    $ snakemake --set-resource-scopes gpu=global disk_mb=local

Resources and Group Jobs
""""""""""""""""""""""""

New to Snakemake 7.11.
When submitting :ref:`group jobs <job_grouping>` to the cluster, Snakemake calculates how many resources to request by first determining which component jobs can be run in parallel, and which must be run in series.
For most resources, such as ``mem_mb`` or ``threads``, a sum will be taken across each parallel layer.
The layer requiring the most resource (i.e. ``max()``) will determine the final amount requested.
The only exception is ``runtime``.
For it, ``max()`` will be used within each layer, then the total amount of time across all layers will be summed.
If resource constraints are provided (via ``--resources`` or ``--cores``) Snakemake will prevent group jobs from requesting more than the constraint.
Jobs that could otherwise be run in parallel will be run in series to prevent the violation of resource constraints.



Preemptible Jobs
""""""""""""""""


You can specify parameters ``preemptible-rules`` and ``preemption-default`` to request a `Google Cloud preemptible virtual machine <https://cloud.google.com/life-sciences/docs/reference/gcloud-examples#using_preemptible_vms>`_ for use with the `Google Life Sciences Executor <https://snakemake.readthedocs.io/en/stable/executing/cloud.html#executing-a-snakemake-workflow-via-google-cloud-life-sciences>`_. There are
several ways to go about doing this. This first example will use preemptible instances for all rules, with 10 repeats (restarts
of the instance if it stops unexpectedly).

.. code-block:: console

    snakemake --preemption-default 10


If your preference is to set a default but then overwrite some rules with a custom value, this is where you can use ``--preemtible-rules``:

.. code-block:: console

    snakemake --preemption-default 10 --preemptible-rules map_reads=3 call_variants=0


The above statement says that we want to use preemtible instances for all steps, defaulting to 10 retries,
but for the steps "map_reads" and "call_variants" we want to apply 3 and 0 retries, respectively. The final
option is to not use preemptible instances by default, but only for a particular rule:


.. code-block:: console

    snakemake --preemptible-rules map_reads=10


Note that this is currently implemented for the Google Life Sciences API.


.. _snakefiles-log:

``log``
^^^^^^

Each rule can specify a log file where information about the execution is written to:

.. code-block:: python

    rule abc:
        input: "input.txt"
        output: "output.txt"
        log: "logs/abc.log"
        shell: "somecommand --log {log} {input} {output}"

Log files can be used as input for other rules, just like any other output file.
However, unlike output files, log files are not deleted upon error.
This is obviously necessary in order to discover causes of errors which might become visible in the log file.

The variable ``log`` can be used inside a shell command to tell the used tool to which file to write the logging information.
The log file has to use the same wildcards as output files, e.g.

.. code-block:: python

    log: "logs/abc.{dataset}.log"


For programs that do not have an explicit ``log`` parameter, you may always use ``2> {log}`` to redirect stderr to a file (here, the ``log`` file) in Linux-based systems.
Note that it is also possible to have multiple named log files, which could be used to capture stdout and stderr:

.. code-block:: python

    rule abc:
        input: "input.txt"
        output: "output.txt"
        log: stdout="logs/foo.stdout", stderr="logs/foo.stderr"
        shell: "somecommand {input} {output} > {log.stdout} 2> {log.stderr}"

``benchmark``
^^^^^^^^^^^^^

Since version 3.1, Snakemake provides support for benchmarking the run times of rules.
This can be used to create complex performance analysis pipelines.
With the `benchmark` keyword, a rule can be declared to store a benchmark of its code into the specified location. E.g. the rule

.. code-block:: python

    rule benchmark_command:
        input:
            "path/to/input.{sample}.txt"
        output:
            "path/to/output.{sample}.txt"
        benchmark:
            "benchmarks/somecommand/{sample}.tsv"
        shell:
            "somecommand {input} {output}"

benchmarks the

* `s`: Wall clock time (in seconds),
* `h:m:s`: Wall clock time (in `hour:minutes:seconds`),
* `max_rss`: Max `RSS <https://en.wikipedia.org/wiki/Resident_set_size>`_ memory usage (in megabytes),
* `max_vms`: Max `VMS <https://en.wikipedia.org/wiki/Virtual_memory>`_ memory usage (in megabytes),
* `max_uss`: Max `USS <https://en.wikipedia.org/wiki/Unique_set_size>`_ memory usage (in megabytes),
* `max_pss`: Max `PSS <https://en.wikipedia.org/wiki/Proportional_set_size>`_ memory usage (in megabytes),
* `io_in`: I/O read (in bytes),
* `io_out`: I/O written (in bytes),
* `mean_load`: CPU load = CPU time (`cpu_usage`) divided by wall clock time (`s`),
* `cpu_time`: CPU time user+system (seconds),

Since version 8.11.0, it is possible to have extra benchmark metrics with the command ``--benchmark-extended``:

* `jobid`: Internal job ID,
* `rule_name`: Rule name,
* `wildcards`: Job wildcards,
* `params`: Job parameters,
* `threads`: Number of threads requested for this job,
* `cpu_usage`: Total CPU load,
* `resources`: Resources requested for this job,
* `input_size_mb`: Size of input files (MiB),

of the command ``somecommand`` for the given output and input files.

For this, the shell or run body of the rule is executed on that data, and all run times are stored into the given benchmark `tsv` file (which will contain a tab-separated table of run times and memory usage in MiB).
Per default, Snakemake executes the job once, generating one run time.
However, the benchmark file can be annotated with the desired number of repeats, e.g.,

.. code-block:: python

    rule benchmark_command:
        input:
            "path/to/input.{sample}.txt"
        output:
            "path/to/output.{sample}.txt"
        benchmark:
            repeat("benchmarks/somecommand/{sample}.tsv", 3)
        shell:
            "somecommand {input} {output}"

will instruct Snakemake to run each job of this rule three times and store all measurements in the benchmark file.
The resulting `tsv` file can be used as input for other rules, just like any other output file.

Since version 8.11.0, it is also possible to have the benchmark metrics in different formats (depending on the extension); currently only the `.jsonl` extension (JSONL format; i.e. one JSON record per line) is supported and all other extensions will be treated as TSV.

.. note::

    Note that benchmarking is only possible in a reliable fashion for subprocesses (thus for tasks run through the ``shell``, ``script``, and ``wrapper`` directive).
    In the ``run`` block, the variable ``bench_record`` is available that you can pass to ``shell()`` as ``bench_record=bench_record``.
    When using ``shell(..., bench_record=bench_record)``, the maximum of all measurements of all ``shell()`` calls will be used but the running time of the rule execution including any Python code.


``messages``
^^^^^^^^^^^^

When executing snakemake, a short summary for each running rule is given to the console. This can be overridden by specifying a message for a rule:


.. code-block:: python

    rule NAME:
        input: "path/to/inputfile", "path/to/other/inputfile"
        output: "path/to/outputfile", "path/to/another/outputfile"
        threads: 8
        message: "Executing somecommand with {threads} threads on the following files {input}."
        shell: "somecommand --threads {threads} {input} {output}"

Note that access to wildcards is also possible via the variable ``wildcards`` (e.g, ``{wildcards.sample}``), which is the same as with shell commands. It is important to have a namespace around wildcards in order to avoid clashes with other variable names.



.. _snakefiles-params:

``params``
^^^^^^^^^^

Sometimes you may want to define certain parameters separately from the rule body. Snakemake provides the ``params`` keyword for this purpose:


.. code-block:: python

    rule:
        input:
            ...
        params:
            threshold=0.4
        output:
            "somedir/{sample}.csv"
        shell:
            "somecommand --threshold {params.threshold} -o {output}"

The ``params`` section is an excellent place to name and assign parameters and variables for your subsequent command.
Similar to ``input``, ``params`` can take functions as well (see :ref:`snakefiles-input_functions`), e.g. you can write

.. code-block:: python

    rule:
        input:
            ...
        params:
            threshold=lambda wildcards: config["thresholds"][wildcards.sample]
        output:
            "somedir/{sample}.csv"
        shell:
            "somecommand --threshold {params.threshold} -o {output}"

Above example mimics a case where one would have to look up the value of some threshold in a config dictionary.
Note that in contrast to the ``input`` directive, functions passed to the
``params`` directive can optionally take more arguments than only ``wildcards``, namely ``input``, ``output``, ``threads``, and ``resources``.
Their order does not matter, apart from the fact that ``wildcards`` has to be the first argument.
This way, params can be used to dynamically adjust those values into whatever format is needed for your command or script.

The ``params`` directive is particularly powerful in combination with Snakemake's :ref:`semantic helper functions <snakefiles-semantic-helpers>`.



.. _snakefiles-paramspace:

Parameter space exploration
"""""""""""""""""""""""""""

The basic Snakemake functionality already provides everything to handle parameter spaces in any way (sub-spacing for certain rules and even depending on wildcard values, the ability to read or generate spaces on the fly or from files via pandas, etc.).
However, it usually would require some boilerplate code for translating a parameter space into wildcard patterns, and translate it back into concrete parameters for scripts and commands.
From Snakemake 5.31 on (inspired by `JUDI <https://pyjudi.readthedocs.io>`_), this is solved via the Paramspace helper, which can be used as follows:

.. code-block:: python

    from snakemake.utils import Paramspace
    import pandas as pd

    # declare a dataframe to be a paramspace
    paramspace = Paramspace(pd.read_csv("params.tsv", sep="\t"))


    rule all:
        input:
            # Aggregate over entire parameter space (or a subset thereof if needed)
            # of course, something like this can happen anywhere in the workflow (not
            # only at the end).
            expand("results/plots/{params}.pdf", params=paramspace.instance_patterns)


    rule simulate:
        output:
            # format a wildcard pattern like "alpha~{alpha}/beta~{beta}/gamma~{gamma}"
            # into a file path, with alpha, beta, gamma being the columns of the data frame
            f"results/simulations/{paramspace.wildcard_pattern}.tsv"
        params:
            # automatically translate the wildcard values into an instance of the param space
            # in the form of a dict (here: {"alpha": ..., "beta": ..., "gamma": ...})
            simulation=paramspace.instance
        script:
            "scripts/simulate.py"


    rule plot:
        input:
            f"results/simulations/{paramspace.wildcard_pattern}.tsv"
        output:
            f"results/plots/{paramspace.wildcard_pattern}.pdf"
        shell:
            "touch {output}"


In above example, **please note** the Python ``f``-string formatting (the ``f`` before the initial quotes) applied to the input and output file strings that contain ``paramspace.wildcard_pattern``.
This means that the file that is registered as input or output file by Snakemake does not contain a wildcard ``{paramspace.wildcard_pattern}``, but instead this item is replaced by a pattern of multiple wildcards derived from the columns of the parameter space dataframe.
This is done by the Python ``f``-string formatting before the string is registered in the rule.
Given that `params.tsv` contains:

.. code-block:: none

    alpha	beta	gamma
    1.0	0.1	0.99
    2.0	0.0	3.9


This workflow will run as follows:

.. code-block:: none

    [Fri Nov 27 20:57:27 2020]
    rule simulate:
        output: results/simulations/alpha~2.0/beta~0.0/gamma~3.9.tsv
        jobid: 4
        wildcards: alpha=2.0, beta=0.0, gamma=3.9

    [Fri Nov 27 20:57:27 2020]
    rule simulate:
        output: results/simulations/alpha~1.0/beta~0.1/gamma~0.99.tsv
        jobid: 2
        wildcards: alpha=1.0, beta=0.1, gamma=0.99

    [Fri Nov 27 20:57:27 2020]
    rule plot:
        input: results/simulations/alpha~2.0/beta~0.0/gamma~3.9.tsv
        output: results/plots/alpha~2.0/beta~0.0/gamma~3.9.pdf
        jobid: 3
        wildcards: alpha=2.0, beta=0.0, gamma=3.9


    [Fri Nov 27 20:57:27 2020]
    rule plot:
        input: results/simulations/alpha~1.0/beta~0.1/gamma~0.99.tsv
        output: results/plots/alpha~1.0/beta~0.1/gamma~0.99.pdf
        jobid: 1
        wildcards: alpha=1.0, beta=0.1, gamma=0.99


    [Fri Nov 27 20:57:27 2020]
    localrule all:
        input: results/plots/alpha~1.0/beta~0.1/gamma~0.99.pdf, results/plots/alpha~2.0/beta~0.0/gamma~3.9.pdf
        jobid: 0


Naturally, it is possible to create sub-spaces from ``Paramspace`` objects, simply by applying all the usual methods and attributes that Pandas data frames provide (e.g. ``.loc[...]``, ``.filter()`` etc.).
Further, the form of the created ``wildcard_pattern`` can be controlled via additional arguments of the ``Paramspace`` `constructor <https://snakemake-api.readthedocs.io/en/latest/api_reference/snakemake_utils.html#snakemake.utils.Paramspace>`_.
In particular, using the argument ``single_wildcard`` the default behavior of encoding each column as a wildcard can be replaced with a single given wildcard name.
This can be handy in case a rule shall serve multiple param spaces with different sets of columns.



.. _snakefiles_retries:

``retries``
^^^^^^^^^^^

Sometimes, rules may be expected to fail occasionally.
For example, this can happen when a rule downloads some online resources.
For such cases, it is possible to defined a number of automatic retries for each job from that particular rule via the ``retries`` directive:

.. code-block:: python

    rule a:
        output:
            "test.txt"
        retries: 3
        shell:
            "curl https://some.unreliable.server/test.txt > {output}"

Often, it is a good idea to combine retry functionality with :ref:`ensure annotations <snakefiles_ensure>`, e.g. for retrying upon invalid checksums or empty files.

Note that it is also possible to define retries globally (via the ``--retries`` command line option, see :ref:`all_options`).
The local definition of the rule thereby overwrites the global definition.

Importantly the ``retries`` directive is meant to be used for defining platform independent behavior (like adding robustness to above download command).
For dealing with unreliable cluster or cloud systems, you should use the ``--retries`` command line option.


.. _snakefiles-job_lifetime_handlers:

``onstart``, ``onsuccess`` and ``onerror``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Sometimes, it is necessary to specify code that shall be executed when the workflow execution is finished (e.g. cleanup, or notification of the user).
With Snakemake 3.2.1, this is possible via the ``onsuccess`` and ``onerror`` keywords:

.. code-block:: python

    onsuccess:
        print("Workflow finished, no error")

    onerror:
        print("An error occurred")
        shell("mail -s "an error occurred" youremail@provider.com < {log}")

The ``onsuccess`` handler is executed if the workflow finished without error. Otherwise, the ``onerror`` handler is executed.
In both handlers, you have access to the variable ``log``, which contains the path to a logfile with the complete Snakemake output.
Snakemake 3.6.0 adds an ``onstart`` handler, that will be executed before the workflow starts.
Note that dry-runs do not trigger any of the handlers.

When you are using :ref:`snakefiles-modules`, only the ``onstart``, ``onsuccess`` and ``onerror`` handlers of the top-level Snakefile are executed. Handlers defined inside module Snakefiles are not triggered automatically.
To access the handlers from a specific module's Snakefile, you can use ``module_name.onstart``, ``module_name.onsuccess`` and ``module_name.onerror``.

.. code-block:: python

    module test1:
        snakefile:
            "module1/Snakefile"

    use rule * from test1 as module1_*

    onstart:
        test1.onstart()
    onsuccess:
        test1.onsuccess()
    onerror:
        test1.onerror()



.. _snakefiles-template-integration:

``template_engine``
^^^^^^^^^^^^^^^^^^^

Sometimes, data analyses entail the dynamic rendering of internal configuration files that are required for certain steps.
From Snakemake 7 on, such template rendering is directly integrated such that it can happen with minimal code and maximum performance.
Consider the following example:

.. code-block:: python

    rule render_jinja2_template:
        input:
            "some-jinja2-template.txt"
        output:
            "results/{sample}.rendered-version.txt"
        params:
            foo=0.1
        template_engine:
            "jinja2"

Here, Snakemake will automatically use the specified template engine `Jinja2 <https://jinja.palletsprojects.com/>`_ to render the template given as input file into the given output file.
The template_engine instruction has to be specified at the end of the rule.
Template rendering rules may only have a single output file.
If the rule needs more than one input file, there has to be one input file called ``template``, pointing to the main template to be used for the rendering:

.. code-block:: python

    rule render_jinja2_template:
        input:
            template="some-jinja2-template.txt",
            other_file="some-other-input-file-used-by-the-template.txt"
        output:
            "results/{sample}.rendered-version.txt"
        params:
            foo=0.1
        template_engine:
            "jinja2"

The template itself has access to ``input``, ``params``, ``wildcards``, and ``config``,
which are the same objects you can use for example in the ``shell`` or ``run`` directive,
and the same objects as can be accessed from ``script`` or ``notebook`` directives (but in the latter two cases they are stored behind the ``snakemake`` object which serves as a dedicated namespace to avoid name clashes).

An example Jinja2 template could look like this::

    This is some text and now we access {{ params.foo }}.

Apart from Jinja2, Snakemake supports `YTE <https://github.com/koesterlab/yte>`_ (YAML template engine), which is particularly designed to support templating of the ubiquitous YAML file format:

.. code-block:: python

    rule render_yte_template:
        input:
            "some-yte-template.yaml"
        output:
            "results/{sample}.rendered-version.yaml"
        params:
            foo=0.1
        template_engine:
            "yte"

Analogously to the jinja2 case YTE has access to ``params``, ``wildcards``, and ``config``:

.. code-block:: yaml

    ?if params.foo < 0.5:
        x:
        - 1
        - 2
        - 3
    ?else:
        y:
        - a
        - b
        - ?config["threshold"]

By default, template rendering rules are executed locally, without submission to cluster or cloud processes (since templating is usually not resource intensive).
However, if a :ref:`storage plugin <storage-support>` is used, a template rule can theoretically leak paths to local copies of the storage files into the rendered template.
This can happen if the template inserts the path of an input file into the rendered output.
Snakemake tries to detect such cases by checking the template output.
To avoid such leaks (only required if your template does something like that with an input file path), you can assign the same :ref:`group <job_grouping>` to your template rule and the consuming rule, and in addition mark the template output as ``temp()``, i.e.:

.. code-block:: python

    rule render_yte_template:
        input:
            "some-yte-template.yaml"
        output:
            temp("results/{sample}.rendered-version.yaml")
        params:
            foo=0.1
        group: "some-group"
        template_engine:
            "yte"

    rule consume_template:
        input:
            "results/{sample}.rendered-version.yaml"
        output:
            "results/some-output.txt"
        group: "some-group"
        shell:
            "sometool {input} {output}"
