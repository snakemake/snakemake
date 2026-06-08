========================
Execution and Scheduling
========================


Local Rules
-----------

When working in a cluster environment, not all rules need to become a job that has to be submitted (e.g. downloading some file, or a target-rule like `all`, see :ref:`snakefiles-targets`).
The keyword `localrules` allows to mark a rule as local, so that it is not submitted to the cluster and instead executed on the host node:

.. code-block:: python

    localrules: all, foo

    rule all:
        input: ...

    rule foo:
        ...

    rule bar:
        ...

Here, only jobs from the rule ``bar`` will be submitted to the cluster, whereas all and foo will be run locally.
Note that you can use the localrules directive **multiple times**. The result will be the union of all declarations.

Alternatively, you can also use the rule directive `localrule`:

.. code-block:: python

    rule all:
        input: ...
        localrule: True

    rule foo:
        ...
        localrule: True

    rule bar:
        ...


        .. _snakefiles-grouping:

Defining groups for execution
-----------------------------

From Snakemake 5.0 on, it is possible to assign rules to groups.
Such groups will be executed together in **cluster** or **cloud mode**, as a so-called **group job**, i.e., all jobs of a particular group will be submitted at once, to the same computing node.
When executing locally, group definitions are ignored.

Groups can be defined via the ``group`` keyword.
This way, queueing and execution time can be saved, in particular if one or several short-running rules are involved.

.. code-block:: python

    samples = [1,2,3,4,5]


    rule all:
        input:
            "test.out"


    rule a:
        output:
            "a/{sample}.out"
        group: "mygroup"
        shell:
            "touch {output}"


    rule b:
        input:
            "a/{sample}.out"
        output:
            "b/{sample}.out"
        group: "mygroup"
        shell:
            "touch {output}"


    rule c:
        input:
            expand("b/{sample}.out", sample=samples)
        output:
            "test.out"
        shell:
            "touch {output}"

Here, jobs from rule ``a`` and ``b`` end up in one group ``mygroup``, whereas jobs from rule ``c`` are executed separately.
Note that Snakemake always determines a **connected subgraph** with the same group id to be a **group job**.
Here, this means that, e.g., the jobs creating ``a/1.out`` and ``b/1.out`` will be in one group, and the jobs creating ``a/2.out`` and ``b/2.out`` will be in a separate group.
However, if we would add ``group: "mygroup"`` to rule ``c``, all jobs would end up in a single group, including the one spawned from rule ``c``, because ``c`` connects all the other jobs.

Alternatively, groups can be defined via the command line interface.
This enables to almost arbitrarily partition the DAG, e.g. in order to save network traffic, see :ref:`here <job_grouping>`.

For execution on the cloud using Google Life Science API and preemptible instances, we expect all rules in the group to be homogeneously set as preemptible instances (e.g., with command-line option ``--preemptible-rules``), such that a preemptible VM is requested for the execution of the group job.

.. _snakefiles_group-local:

Group-local jobs
~~~~~~~~~~~~~~~~

From Snakemake 7.0 on, it is further possible to ensure that jobs from a certain rule are executed separately within each :ref:`job group <job_grouping>`.
For this purpose we use :ref:`input functions <snakefiles-input_functions>`, which, in addition to the ``wildcards`` argument can expect a ``groupid`` argument.
In such a case, Snakemake passes the ID of the corresponding group job to the input function.
Consider the following example

.. code-block:: python

    rule all:
        input:
            expand("bar{i}.txt", i=range(3))


    rule grouplocal:
        output:
            "foo.{groupid}.txt"
        group:
            "foo"
        shell:
            "echo test > {output}"


    def get_input(wildcards, groupid):
        return f"foo.{groupid}.txt"


    rule consumer:
        input:
            get_input
        output:
            "bar{i}.txt"
        group:
            "foo"
        shell:
            "cp {input} {output}"

Here, the value of ``groupid`` that is passed by Snakemake to the input function is a `UUID <https://en.wikipedia.org/wiki/Universally_unique_identifier>`_ that uniquely identifies the group job in which each instance of the rule ``consumer`` is contained.
In the input function ``get_input`` we use this ID to request the desired input file from the rule ``grouplocal``.
Since the value of the corresponding wildcard ``groupid`` is now always a group specific unique ID, it is ensured that the rule ``grouplocal`` will run for every group job spawned from the group ``foo`` (remember that group jobs by default only span one connected component, and that this can be configured via the command line, see :ref:`job_grouping`).
Of course, above example would also work if the groups are not specified via the rule definition but entirely via the :ref:`command line <job_grouping>`.



.. _snakefiles-service-rules:

Service rules/jobs
------------------

From Snakemake 7.0 on, it is possible to define so-called service rules.
Jobs spawned from such rules provide at least one special output file that is marked as ``service``, which means that it is considered to provide a resource that shall be kept available until all consuming jobs are finished.
This can for example be the socket of a database, a shared memory device, a ramdisk, and so on.
It can even just be a dummy file, and access to the service might happen via a different channel (e.g. a local http port).
Service jobs are expected to not exit after creating that resource, but instead wait until Snakemake terminates them (e.g. via SIGTERM on Unixoid systems).

Consider the following example:

.. code-block:: python

    rule the_service:
        output:
            service("foo.socket")
        shell:
            # here we simulate some kind of server process that provides data via a socket
            "ln -s /dev/random {output}; sleep 10000"


    rule consumer1:
        input:
            "foo.socket"
        output:
            "test.txt"
        shell:
            "head -n1 {input} > {output}"


    rule consumer2:
        input:
            "foo.socket"
        output:
            "test2.txt"
        shell:
            "head -n1 {input} > {output}"

Snakemake will schedule the service with all consumers to the same physical node (in the future we might provide further controls and other modes of operation).
Once all consumer jobs are finished, the service job will be terminated automatically by Snakemake, and the service output will be removed.

Group-local service jobs
^^^^^^^^^^^^^^^^^^^^^^^^

Since Snakemake supports arbitrary partitioning of the DAG into so-called :ref:`job groups <job_grouping>`, one should consider what this implies for service jobs when running a workflow in a cluster of cloud context:
since each group job spans at least one connected component (see :ref:`job groups <job_grouping>` and `the Snakemake paper <https://doi.org/10.12688/f1000research.29032.2>`), this means that the service job will automatically connect all consumers into one big group.
This can be undesired, because depending on the number of consumers that group job can become too big for efficient execution on the underlying architecture.
In case of local execution, this is not a problem because here DAG partitioning has no effect.

However, to make a workflow portable across different backends, this behavior should always be considered.
In order to circumvent it, it is possible to model service jobs as group-local, i.e. ensuring that each group job gets its own instance of the service rule.
This works by combining the service job pattern from above with the :ref:`group-local pattern <snakefiles_group-local>` as follows:

.. code-block:: python

    rule the_service:
        output:
            service("foo.{groupid}.socket")
        shell:
            # here we simulate some kind of server process that provides data via a socket
            "ln -s /dev/random {output}; sleep 10000"


    def get_socket(wildcards, groupid):
        return f"foo.{groupid}.socket"


    rule consumer1:
        input:
            get_socket
        output:
            "test.txt"
        shell:
            "head -n1 {input} > {output}"


    rule consumer2:
        input:
            get_socket
        output:
            "test2.txt"
        shell:
            "head -n1 {input} > {output}"



Priorities
----------

Snakemake allows for rules that specify numeric and/or callable priorities:


.. code-block:: python

    rule:
        input: ...
        output: ...
        priority: 50
        shell: ...

Per default, each rule has a priority of 0. Any rule that specifies a higher priority, will be preferred by the scheduler over all rules that are ready to execute at the same time without having at least the same priority.

Priority may also be specified with a callable. The callable receives ``wildcards``
as its first positional argument, and may optionally accept ``input``, ``attempt``,
and ``rulename`` keyword arguments (similar to param functions). It has to return the priority as an integer or float (will be rounded):

.. code-block:: python

    rule sort:
        input: "{dataset}.txt"
        output: "{dataset}.sorted"
        priority: lambda wildcards, input: input.size_mb
        shell: "sort {input} > {output}"

This allows the scheduler to dynamically prioritise jobs based on, e.g., input
file size so that larger jobs start first.

Furthermore, the ``--prioritize`` or ``-P`` command line flag allows to specify files (or rules) that shall be created with highest priority during the workflow execution. This means that the scheduler will assign the specified target and all its dependencies highest priority, such that the target is finished as soon as possible.
The ``--dry-run`` (equivalently ``--dryrun``) or ``-n`` option allows you to see the scheduling plan including the assigned priorities.
