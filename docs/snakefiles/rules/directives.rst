
.. _snakefiles-directives:

Directives
----------

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




.. _snakefiles-log:

``log``
^^^^^^^

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



.. _snakefiles-threads:

``threads``
^^^^^^^^^^^

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

