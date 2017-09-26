.. _executable:

===================
Executing Snakemake
===================

This part of the documentation describes the ``snakemake`` executable.  Snakemake
is primarily a command-line tool, so the ``snakemake`` executable is the primary way
to execute, debug, and visualize workflows.

.. user_manual-snakemake_options:

-----------------------------
Useful Command Line Arguments
-----------------------------

If called without parameters, i.e.

.. code-block:: console

    $ snakemake

Snakemake tries to execute the workflow specified in a file called ``Snakefile`` in the same directory (instead, the Snakefile can be given via the parameter ``-s``).

By issuing

.. code-block:: console

    $ snakemake -n

a dry-run can be performed.
This is useful to test if the workflow is defined properly and to estimate the amount of needed computation.
Further, the reason for each rule execution can be printed via


.. code-block:: console

    $ snakemake -n -r

Importantly, Snakemake can automatically determine which parts of the workflow can be run in parallel.
By specifying the number of available cores, i.e.

.. code-block:: console

    $ snakemake -j 4

one can tell Snakemake to use up to 4 cores and solve a binary knapsack problem to optimize the scheduling of jobs.
If the number is omitted (i.e., only ``-j`` is given), the number of used cores is determined as the number of available CPU cores in the machine.


-------------
Cloud Support
-------------

Snakemake 4.0 and later supports experimental execution in the cloud via Kubernetes.
This is independent of the cloud provider, but we provide the setup steps for GCE below.

Google cloud engine
~~~~~~~~~~~~~~~~~~~

First, install the `Google Cloud SDK <https://cloud.google.com/sdk/docs/quickstarts>`_.
Then, run

.. code-block:: console

    $ gcloud init

to setup your access.
Then, you can create a new kubernetes cluster via

.. code-block:: console

    $ gcloud container clusters create $CLUSTER_NAME --num-nodes=$NODES --scopes storage-rw

with ``$CLUSTER_NAME`` being the cluster name and ``$NODES`` being the number of cluster
nodes. If you intent to use google storage, make sure that `--scopes storage-rw` is set.
This enables Snakemake to write to the google storage from within the cloud nodes.
Next, you configure Kubernetes to use the new cluster via

.. code-block:: console

    $ gcloud container clusters get-credentials $CLUSTER_NAME


Now, Snakemake is ready to use your cluster.

**Important:** After finishing your work, do not forget to delete the cluster with

.. code-block:: console

    $ gcloud container clusters delete $CLUSTER_NAME

in order to avoid unnecessary charges.


.. _kubernetes:

Executing a Snakemake workflow via kubernetes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Assuming that kubernetes has been properly configured (see above), you can
execute a workflow via:

.. code-block:: console

    snakemake --kubernetes --use-conda --default-remote-provider $REMOTE --default-remote-prefix $PREFIX

In this mode, Snakemake will assume all input and output files to be stored in a given
remote location, configured by setting ``$REMOTE`` to your provider of choice
(e.g. ``GS`` for Google cloud storage or ``S3`` for Amazon S3) and ``$PREFIX``
to a bucket name or subfolder within that remote storage.
After successful execution, you find your results in the specified remote storage.
Of course, if any input or output already defines a different remote location, the latter will be used instead.
Importantly, this means that Snakemake does **not** require a shared network
filesystem to work in the cloud.

It is further possible to forward arbitrary environment variables to the kubernetes
jobs via the flag ``--kubernetes-env`` (see ``snakemake --help``).

When executing, Snakemake will make use of the defined resources and threads
to schedule jobs to the correct nodes. In particular, it will forward memory requirements
defined as `mem_mb` to kubernetes. Further, it will propagate the number of threads
a job intends to use, such that kubernetes can allocate it to the correct cloud
computing node.


-----------------
Cluster Execution
-----------------


Snakemake can make use of cluster engines that support shell scripts and have access to a common filesystem, (e.g. the Sun Grid Engine).
In this case, Snakemake simply needs to be given a submit command that accepts a shell script as first positional argument:

.. code-block:: console

    $ snakemake --cluster qsub -j 32


Here, ``-j`` denotes the number of jobs submitted being submitted to the cluster at the same time (here 32).
The cluster command can be decorated with job specific information, e.g.

.. code-block:: console

    $ snakemake --cluster "qsub {threads}"

Thereby, all keywords of a rule are allowed (e.g. params, input, output, threads, priority, ...).
For example, you could encode the expected running time into params:

.. code-block:: python

    rule:
        input:  ...
        output: ...
        params: runtime="4h"
        shell: ...

and forward it to the cluster scheduler:

.. code-block:: console

    $ snakemake --cluster "qsub --runtime {params.runtime}"

If your cluster system supports `DRMAA <http://www.drmaa.org/>`_, Snakemake can make use of that to increase the control over jobs.
E.g. jobs can be cancelled upon pressing ``Ctrl+C``, which is not possible with the generic ``--cluster`` support.
With DRMAA, no ``qsub`` command needs to be provided, but system specific arguments can still be given as a string, e.g.

.. code-block:: console

    $ snakemake --drmaa " -q username" -j 32

Note that the string has to contain a leading whitespace.
Else, the arguments will be interpreted as part of the normal Snakemake arguments, and execution will fail.

Job Properties
..............

When executing a workflow on a cluster using the ``--cluster`` parameter (see below), Snakemake creates a job script for each job to execute. This script is then invoked using the provided cluster submission command (e.g. ``qsub``). Sometimes you want to provide a custom wrapper for the cluster submission command that decides about additional parameters. As this might be based on properties of the job, Snakemake stores the job properties (e.g. rule name, threads, input files, params etc.) as JSON inside the job script. For convenience, there exists a parser function `snakemake.utils.read_job_properties` that can be used to access the properties. The following shows an example job submission wrapper:

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


--------
Profiles
--------

Adapting Snakemake to a particular environment can entail many flags and options.
Therefore, since Snakemake 4.1, it is possible to specify a configuration profile
to be used to obtain default options:

.. code-block:: console

   $ snakemake --profile myprofile

Here, a folder ``myprofile`` is searched in per-user and global configuration directories (on Linux, this will be ``$HOME/.config/snakemake`` and ``/etc/xdg/snakemake``, you can find the answer for your system via ``snakemake --help``).
Alternatively, an absolute or relative path to the folder can be given.
The profile folder is expected to contain a file ``config.yaml`` that defines default values for the Snakemake command line arguments.
For example, the file

.. code-block:: yaml

    cluster: qsub
    jobs: 100

would setup Snakemake to always submit to the cluster via the ``qsub`` command, and never use more than 100 parallel jobs in total.
Under https://github.com/snakemake-profiles/doc, you can find publicly available profiles.
Feel free to contribute your own.

The profile folder can additionally contain auxilliary files, e.g., jobscripts, or any kind of wrappers.
See https://github.com/snakemake-profiles/doc for examples.

.. _getting_started-visualization:

-------------
Visualization
-------------

To visualize the workflow, one can use the option ``--dag``.
This creates a representation of the DAG in the graphviz dot language which has to be postprocessed by the graphviz tool ``dot``.
E.g. to visualize the DAG that would be executed, you can issue:

.. code-block:: console

    $ snakemake --dag | dot | display

For saving this to a file, you can specify the desired format:

.. code-block:: console

    $ snakemake --dag | dot -Tpdf > dag.pdf

To visualize the whole DAG regardless of the eventual presence of files, the ``forceall`` option can be used:

.. code-block:: console

    $ snakemake --forceall --dag | dot -Tpdf > dag.pdf

Of course the visual appearance can be modified by providing further command line arguments to ``dot``.


.. _all_options:

-----------
All Options
-----------

.. argparse::
   :module: snakemake
   :func: get_argument_parser
   :prog: snakemake

   All command line options can be printed by calling ``snakemake -h``.

.. _getting_started-bash_completion:

---------------
Bash Completion
---------------

Snakemake supports bash completion for filenames, rulenames and arguments.
To enable it globally, just append

.. code-block:: bash

    `snakemake --bash-completion`

including the accents to your ``.bashrc``.
This only works if the ``snakemake`` command is in your path.
