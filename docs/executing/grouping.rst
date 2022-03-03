.. _job_grouping:

============
Job Grouping
============

The graph of jobs that Snakemake determines before execution can be partitioned into groups.
Such groups will be executed together in **cluster** or **cloud mode**, as a so-called **group job**, i.e., all jobs of a particular group will be submitted at once, to the same computing node.
When executing locally, group definitions are ignored.

Groups can be defined along with the workflow definition via the ``group`` keyword, see :ref:`snakefiles-grouping`.
This way, queueing and execution time can be saved, in particular by attaching short-running downstream jobs to long running upstream jobs.

However, often the benefit grouping should be heavily dependent on the specifics of the underlying computing platform.
Hence, it is possible to assign groups via the command line.
For example, with

.. code-block:: bash

    snakemake --groups somerule=group0 someotherrule=group0

we assign the two rules ``somerule`` and ``someotherrule`` to the same group ``group0``.
by default, groups do not span disconnected parts of the DAG.
Here, this means that by default only jobs of ``somerule`` and ``someotherrule`` end in the same group that are directly connected.
It is however possible to configure the number of connected DAG components that are spanned by a group via the flag ``--group-components``.
This way, it is e.g. possible to define batches of jobs of the same kind that shall be executed within one group, e.g.


.. code-block:: bash

    snakemake --groups somerule=group0 --group-components group0=5

means that given that there exist ``n`` jobs spawned from rule ``somerule``, Snakemake will create ``n / 5`` groups which each execute 5 jobs of ``somerule`` together.
For example, with 10 jobs from ``somerule`` you would end up with 2 groups of 5 jobs that are submitted as one piece each.
