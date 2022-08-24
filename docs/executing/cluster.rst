.. _cluster:

=================
Cluster Execution
=================


Snakemake can make use of cluster engines that support shell scripts and have access to a common filesystem, (e.g. the Sun Grid Engine).
In this case, Snakemake simply needs to be given a submit command that accepts a shell script as first positional argument:

.. code-block:: console

    $ snakemake --cluster qsub --jobs 32


Here, ``--jobs`` denotes the number of jobs submitted to the cluster at the same time (here 32).
The cluster command can be decorated with job specific information, e.g.

.. sidebar:: Note

  Consider to :ref:`group jobs <snakefiles-grouping>` in order to minimize overhead, in particular for short-running jobs.


.. code-block:: console

    $ snakemake --cluster "qsub {threads}"

Thereby, all keywords of a rule are allowed (e.g. rulename, params, input, output, threads, priority, resources, ...).
For example, you could encode the expected running time in minutes into a :ref:`resource <snakefiles-resources>` ``runtime_min``:

.. code-block:: python

    rule:
        input:  
            ...
        output:
            ...
        resources: 
            runtime_min=240
        shell:
            ...

and forward it to the cluster scheduler:

.. code-block:: console

    $ snakemake --cluster "qsub --runtime {resources.runtime}"

In order to avoid specifying ``runtime_min`` for each rule, you can make use of the ``--default-resources`` flag, see ``snakemake --help``.

If your cluster system supports `DRMAA <https://www.drmaa.org/>`_, Snakemake can make use of that to control jobs.
With DRMAA, no ``qsub`` command needs to be provided, but system specific arguments can still be given as a string, e.g.

.. code-block:: console

    $ snakemake --drmaa " -q username" -j 32

Note that the string has to contain a leading whitespace.
Else, the arguments will be interpreted as part of the normal Snakemake arguments, and execution will fail.

Adapting to a specific cluster can involve quite a lot of options. It is therefore a good idea to setup a :ref:`a profile <profiles>`.

--------------
Executing on SLURM clusters
--------------

`SLURM <https://slurm.schedmd.com/documentation.html>(the Simple Linux Utility for Resource Management)` is one a widely used batch system for
performance compute clusters. In order to use Snakemake with slurm, simply append `--slurm` to your command line.

- Specifying `Account` and `Partition`
  
  Most SLURM clusters have two mandatory resource indicators for accounting and scheduling, `Account` and `Partition`, respectivily.
  These resources are usually omitted from Snakemake workflows. 

  To specify them issue a default:

  .. code-block:: console

    $ snakemake --slurm --default-resources account=<your SLURM account> 


- General Resource Specifications:

  A workflow rule may support a number of :ref:`resource <snakefiles-resources>` specification. For a SLURM cluster, 
  a mapping between Snakemake and SLURM needs to be performed.

  We can use the following specifications, unique per rule:

  +-----------------+-----------------------+------------------------------------------------------------------+
  | SLURM Resource  | Snakemake Equivalent  | Background Information                                           |
  +=================+=======================+==================================================================+
  | `-p`/`--partition`   | `partition`      | the partition a rule/job is to use                               |
  +-----------------+-----------------------+------------------------------------------------------------------+
  | `-t`/`--time`   | `walltime_minutes`    | the walltime per job in minutes                                  |
  +-----------------+-----------------------+------------------------------------------------------------------+
  | `-C`/`constraint`| `constraint`         | may hold features on some clusters                               |
  +-----------------+-----------------------+------------------------------------------------------------------+
  | `--mem[unit]`   |  `mem_mb`             | memory in `unit` a cluster node must provide                     |
  +-----------------+-----------------------+------------------------------------------------------------------+
  | `--mem-per-cpu` |  `mem_mb_per_cpu`     | memory per reserved CPU                                          |
  +-----------------+-----------------------+------------------------------------------------------------------+
  |  `-n`/`ntasks`  |  `ntasks`             | number of concurrent tasks / ranks                               |
  +-----------------+-----------------------+------------------------------------------------------------------+
  | `-c`/`--cpus-per-taks` | `cpus_per_task`| number of cpus per task                                          |
  +-----------------+-----------------------+------------------------------------------------------------------+
  | `-N`/`--nodes`  | `nodes`               | number of nodes                                                  |
  +-----------------+-----------------------+------------------------------------------------------------------+
  
  Each of these parameters can be part of a rule, e.g.:

  .. code-block:: python
    rule:
      input: ...
      output: ...
      resources:
          partition: <partition name>
          walltime_minutes: <some number>

  Please note: as `--mem` and `--mem-per-cpu` are mutually exclusive on SLURM clusters, there corresponding resource flags `mem_mb` and `mem_mb_per_cpu` are mutually exclusive, too.
  You can only reserve memory a compute node has to provide or the memory required per CPU (SLURM does not make any distintion between real CPU cores and those provided by
  hyperthreads). SLURM will try to sastify a combination of `mem_mb_per_cpu` and `cpus_per_task` and `nodes`, if `nodes` is not given.

- Ordinary SMP jobs
  
  Most jobs will be carried out by programs which are either single core scripts or threaded programs, hence SMP (:ref: shared memory programs<https://en.wikipedia.org/wiki/Shared_memory>)
  in nature. Here, at most `mem_mb_per_cpu` and `cpus_per_task` need to be specified, e.g.:

  .. code-block:: python
    rule:
      input: ...
      output: ...
      resources:
        mem_mb_per_cpu: 1800
        cpus_per_task: 8

  This will give your application 8*1800M, hence 14.4 GB. Please consult your site's documentation for your batch system default settings and memory distribution.

- Group Jobs

  :ref: group jobs<snakefiles-grouping> may not oversubscribe one compute node. Best simply ask for one compute node with sufficient memory, e.g.:
  
  .. code-block:: python
    rule:
      input: ...
      output: ...
      resources:
        mem_mb: 57000 # or more, consult your site's doctumentation
        nodes: 1

- MPI jobs
  
  Highly parallel programs often use the MPI (:ref: message passing interface<https://en.wikipedia.org/wiki/Message_Passing_Interface>) to enable a programm to span
  and work across an invidual compute node's boundary. Und SLURM `srun` is the MPI-starter. Snakemake supports using MPI-programs like this:

  .. code-block:: python
    rule:
      input: ...
      output: ...
      resources:
        mpi: srun

  Here, `srun` may be followed by other instructions to get your communication or MPI-topology right, e.g.:

  .. code-block:: python
      resources:
        mpi: srun --cpu-bind=ranks

  or

  .. code-block:: python
      resources:
        mpi: srun --hint=nomultithread --mpi=pmix

  In any case Snakemake will start your MPI-program with this entire line after the `mpi`-directive. You may also choose a different MPI-starter, e.g. `mpirun`, instead.



--------------
Job Properties
--------------

When executing a workflow on a cluster using the ``--cluster`` parameter (see below), Snakemake creates a job script for each job to execute. This script is then invoked using the provided cluster submission command (e.g. ``qsub``). Sometimes you want to provide a custom wrapper for the cluster submission command that decides about additional parameters. As this might be based on properties of the job, Snakemake stores the job properties (e.g. name, rulename, threads, input, output, params etc.) as JSON inside the job script (for group jobs, the rulename will be "GROUP", otherwise it will be the same as the job name). For convenience, there exists a parser function `snakemake.utils.read_job_properties` that can be used to access the properties. The following shows an example job submission wrapper:

.. code-block:: python

    #!python

    #!/usr/bin/env python3
    import os
    import sys

    from snakemake.utils import read_job_properties

    jobscript = sys.argv[1]
    job_properties = read_job_properties(jobscript)

    # do something useful with the threads
    threads = job_properties[threads]

    # access property defined in the cluster configuration file (Snakemake >=3.6.0)
    job_properties["cluster"]["time"]

    os.system("qsub -t {threads} {script}".format(threads=threads, script=jobscript))
