.. _job_grouping:

============
Job Grouping
============

The graph of jobs that Snakemake determines before execution can be partitioned into groups.
Such groups will be executed together in **cluster** or **cloud mode**, as a so-called **group job**, i.e., all jobs of a particular group will be submitted at once, to the same computing node.
When executing locally, group definitions are ignored.

Groups can be defined along with the workflow definition via the ``group`` keyword, see :ref:`snakefiles-grouping`.
This way, queueing and execution time can be saved, in particular by attaching short-running downstream jobs to long running upstream jobs.

From Snakemake 7.11 on, Snakemake will request resources for groups by summing across jobs that can be run in parallel, and taking the max of jobs run in series.
The only exception is ``runtime``, where the max will be taken over parallel jobs, and the sum over series.
If resource constraints are provided (via ``--resources`` or ``--cores``), parallel job layers that exceed the constraints will be stacked in series.
For example, if 6 instances of ``somerule`` are being run, each instance requires ``1000MB`` of memory and ``30 min`` runtime, and only ``3000MB`` are available, Snakemake will request ``3000MB`` and ``60 min`` runtime, enough to run 3 instances of ``somerule``, then another 3 instances of ``somerule`` in series.

Often, the ideal group will be dependent on the specifics of the underlying computing platform.
Hence, it is possible to assign groups via the command line.
For example, with

.. code-block:: bash

    snakemake --groups somerule=group0 someotherrule=group0

we assign the two rules ``somerule`` and ``someotherrule`` to the same group ``group0``.

By default, groups do not span disconnected parts of the DAG.
This means that, for example, jobs of ``somerule`` and ``someotherrule`` only end in the same group if they are directly connected.
It is, however, possible to configure the number of connected DAG components that are spanned by a group via the flag ``--group-components``.
This makes it possible to define batches of jobs of the same kind that shall be executed within one group. For instance:


.. code-block:: bash

    snakemake --groups somerule=group0 --group-components group0=5

means that given ``n`` jobs spawned from rule ``somerule``, Snakemake will create ``n / 5`` groups which each execute 5 jobs of ``somerule`` together.
For example, with 10 jobs from ``somerule`` you would end up with 2 groups of 5 jobs that are submitted as one piece each.

Furthermore, it is possible to use wildcards in group names.
This way, you can e.g. have a group per sample, e.g.:

.. code-block:: bash

    snakemake --groups somerule=group_{sample} --group-components group_{sample}=5
