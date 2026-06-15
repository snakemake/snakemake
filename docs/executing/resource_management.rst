===================
Resource Management
===================


In general, resources are just names to the Snakemake scheduler, i.e., Snakemake does not check on the resource consumption of jobs in real time.
Instead, resources are used to determine which jobs can be executed at the same time without exceeding the limits specified at the command line.
Apart from making Snakemake aware of hybrid-computing architectures (e.g. with a limited number of additional devices like GPUs) this allows us to control scheduling in various ways, e.g. to limit IO-heavy jobs by assigning an artificial IO-resource to them and limiting it via the ``--resources`` flag.
If no limits are given, the resources are ignored in local execution.

Resources can have any arbitrary name, and must be assigned ``int`` or ``str`` values.
In case of ``None``, the resource is considered to be unset (i.e. ignored) in the rule.

.. _snakefiles-dynamic-resources:

-----------------
Dynamic Resources
-----------------

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

------------------
Standard Resources
------------------

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

-----------------
Default Resources
-----------------

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

------------------------------
Resources and Remote Execution
------------------------------

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

------------------------
Resources and Group Jobs
------------------------

New to Snakemake 7.11.
When submitting :ref:`group jobs <job_grouping>` to the cluster, Snakemake calculates how many resources to request by first determining which component jobs can be run in parallel, and which must be run in series.
For most resources, such as ``mem_mb`` or ``threads``, a sum will be taken across each parallel layer.
The layer requiring the most resource (i.e. ``max()``) will determine the final amount requested.
The only exception is ``runtime``.
For it, ``max()`` will be used within each layer, then the total amount of time across all layers will be summed.
If resource constraints are provided (via ``--resources`` or ``--cores``) Snakemake will prevent group jobs from requesting more than the constraint.
Jobs that could otherwise be run in parallel will be run in series to prevent the violation of resource constraints.


----------------
Preemptible Jobs
----------------


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

