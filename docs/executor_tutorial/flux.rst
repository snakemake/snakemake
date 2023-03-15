
.. _tutorial-flux:

Flux Tutorials
==============

These tutorials will cover using the Flux and Flux Operator executors.

.. _Snakemake: http://snakemake.readthedocs.io
.. _Snakemake Remotes: https://snakemake.readthedocs.io/en/stable/snakefiles/remote_files.html
.. _Python: https://www.python.org/

`Flux-framework <https://flux-framework.org/>`_ is a flexible resource scheduler that can work
on both high performance computing systems and cloud (e.g., Kubernetes).
Since it is more modern (e.g., has an official Python API) we define it under a cloud resource. 
This also means that here we provide two tutorials, and two modes of execution:

 - Local Flux Executor: runs on a local cluster with a Flux Framework instance 
 - Flux Operator: spins up "MiniClusters" on demand to achieve the same


.. _tutorial-flux-executor:

Flux Executor
-------------

The Flux Executor will issue commands to a Flux Scheduler that is locally accessible
via the Flux Python bindings.

Setup
:::::

To go through this tutorial, you need the following software installed:

- Docker

For this example, we will show you how to set up a "single node" local Flux container to interact with snakemake. You can use the `Dockerfile in examples/flux <https://github.com/snakemake/snakemake/blob/main/examples/flux/Dockerfile>`_ that will provide a container with Flux and snakemake        
Note that we install from source and bind to ``/home/fluxuser/snakemake`` with the intention of being able to develop (if desired).
First, build the container:

.. code-block:: console

    $ docker build -t flux-snake .

And then you can run the container with or without any such bind:

.. code-block:: console

    $ docker run -it --rm flux-snake 

Once you shelled into the container, you can view and start a Flux instance:

.. code-block:: console

    $ flux getattr size 
    $ flux start --test-size=4
    $ flux getattr size 

And see resources available:

.. code-block:: console

    $ flux resource status
    STATUS NNODES NODELIST
     avail      4 5a74dc238d[98,98,98,98]


Resources
:::::::::

Flux currently has support for ``runtime``, which should be set to a number
(seconds) and defaults to 0, meaning unlimited runtime. Flux currently does not
support ``mem_mb`` or ``disk_mb``.

         

Run Snakemake
:::::::::::::

Now let's run Snakemake with the Flux executor. There is an example ``Snakefile``
in the flux examples folder that will show running a "Hello World!" example,
and this file should be in your fluxuser home:


.. code-block:: console

    $ ls
    Snakefile  snakemake

Here is how to run the workflow:

.. code:: console

    $ snakemake --flux --jobs=1

The flags above refer to:

 - ``--flux``: tell Snakemake to use the flux executor

Once you submit the job, you'll immediately see the familiar Snakemake console output.
The jobs happen very quickly, but the default wait time between checks is 10 seconds
so it will take a bit longer.

.. code:: console

    Building DAG of jobs...
    Using shell: /usr/bin/bash
    Provided cores: 1 (use --cores to define parallelism)
    Rules claiming more threads will be scaled down.
    Job stats:
    job                         count    min threads    max threads
    ------------------------  -------  -------------  -------------
    all                             1              1              1
    multilingual_hello_world        2              1              1
    total                           3              1              1

    Select jobs to execute...

    [Fri Aug 12 21:09:32 2022]
    rule multilingual_hello_world:
        output: hello/world.txt
        jobid: 1
        reason: Missing output files: hello/world.txt
        wildcards: greeting=hello
        resources: tmpdir=/tmp

    Checking status for job ƒ3sWJLhD
    [Fri Aug 12 21:09:42 2022]
    Finished job 1.
    1 of 3 steps (33%) done
    Select jobs to execute...

    [Fri Aug 12 21:09:42 2022]
    rule multilingual_hello_world:
        output: hola/world.txt
        jobid: 2
        reason: Missing output files: hola/world.txt
        wildcards: greeting=hola
        resources: tmpdir=/tmp

    Checking status for job ƒ8JAY1Kd
    [Fri Aug 12 21:09:52 2022]
    Finished job 2.
    2 of 3 steps (67%) done
    Select jobs to execute...

    [Fri Aug 12 21:09:52 2022]
    localrule all:
        input: hello/world.txt, hola/world.txt
        jobid: 0
        reason: Input files updated by another job: hola/world.txt, hello/world.txt
        resources: tmpdir=/tmp

    [Fri Aug 12 21:09:52 2022]
    Finished job 0.
    3 of 3 steps (100%) done
    Complete log: .snakemake/log/2022-08-12T210932.564786.snakemake.log

At this point you can inspect the local directory to see your job output!

.. code:: console

    $ ls
    Snakefile  hello  hola
    $ cat hello/world.txt 
    hello, World!


Flux Without Shared Filesystem
::::::::::::::::::::::::::::::

By default, the Flux executor assumes a shared filesystem. If this isn't the case, you can add
the ``--no-shared-fs`` flag, which will tell Snakemake that Flux is running without a shared filesystem.

.. code:: console

    $ snakemake --flux --jobs=1 --no-shared-fs


See the `flux documentation <https://flux-framework.readthedocs.io/en/latest/quickstart.html#docker-recommended-for-quick-single-node-deployments>`_
for more detail. For now, let's try interacting with flux via snakemake via the `Flux Python Bindings <https://flux-framework.readthedocs.io/projects/flux-workflow-examples/en/latest/job-submit-api/README.html>`_.

The code for this example is provided in  (`examples/flux <https://github.com/snakemake/snakemake/tree/main/examples/flux>`_)


.. _tutorial-flux-operator-executor:

Flux Operator Executor
----------------------

The Flux Operator Executor will issue commands to a Kubernetes cluster to run your jobs on a set of
networked pods with a Flux instance running called a 
`MiniCluster <https://flux-framework.org/flux-operator/getting_started/custom-resource-definition.html>`_.

Setup
:::::

To go through this tutorial, you need the following additional software installed or accessible

- Kubernetes (e.g., MiniKube or similar)
- kubectl 

For this tutorial we will show you how to launch MiniKube and create a MiniCluster to run your jobs. We will
do this in two parts - the first assuming that all job steps require the same resources (and thus can run on the
same cluster) and the second assuming the steps need different resources, and we will create more than one
MiniCluster.  First, bring up MiniKube and create a namespace for your jobs:

.. code-block:: console

    $ minikube start
    $ kubectl create namespace flux-operator

Next, prepare the Snakemake tutorial data in a temporary directory, ``/tmp/workflow``.

.. code-block:: console

    $ git clone --depth 1 https://github.com/snakemake/snakemake-tutorial-data /tmp/workflow
    $ mkdir -p /tmp/workflow/scripts
    $ wget -O /tmp/workflow/scripts/plot-quals.py https://raw.githubusercontent.com/rse-ops/flux-hpc/main/snakemake/atacseq/scripts/plot-quals.py


The Snakefile can be found in the ``./examples/flux/operator`` directory of Snakemake:

.. code-block:: console

    $ cp ./examples/flux/operator/Snakefile /tmp/workflow/Snakefile

The main difference is that it has a container defined for each step. Since we are using MiniKube
(which isn't great at pulling containers like a production cluster) let's pull the container first:

.. code-block:: console

    $ minikube ssh docker pull ghcr.io/rse-ops/atacseq-vanilla:app-latest

And finally, in a separate terminal, make sure your host ``/tmp/workflow`` is bound to the same
path in the MiniKube virtual machine:

.. code-block:: console

    $ minikube ssh -- mkdir -p /tmp/workflow
    $ minikube mount /tmp/workflow:/tmp/workflow 

You'll need to install the Flux Operator! This is the easiest way:

.. code-block:: console

    $ wget https://raw.githubusercontent.com/flux-framework/flux-operator/main/examples/dist/flux-operator.yaml
    $ kubectl apply -f flux-operator.yaml 

Run the Workflow 
::::::::::::::::

Finally, run the workflow from this directory, and ask for the flux Operator.
Note that this currently depends on having the same input files locally and 
in the cluster (created as a volume), and this is because we are running the demo 
locally. For a production run, we'd want to mount the volume from cloud storage,
and then run the command from there.

.. code-block:: console

    $ snakemake --cores 1 --jobs 1 --flux-operator

And you'll see the jobs run! When it's done (and you see outputs) try deleting everything,
and then running again and allowing for more than one job to be run at once.

.. code-block:: console

    $ rm -rf calls/ mapped_reads/ sorted_reads/ plots/
    $ snakemake --cores 1 --jobs 2 --flux-operator

You'll notice the workflow moving faster, and this is because we have submit more
than one job at once!

How does it work?
:::::::::::::::::

MiniClusters are unique based on the container base needed, and the number of nodes for a specific level,
where the level include a set of jobs that can be run in parallel. Since we are always paying for the
Kubernetes cluster, no matter what, we will do best to optimize its usage, meaning running as many
jobs in parallel at once as we can. For this reason, we use a simple algorithm that calculates
the maximum number of nodes that might be used at a single level (up to the max size of your Kubernetes cluster)
and use that paired with the container image to define a unique ID for a cluster. This means that when you are 
using the Flux Operator executor you should:

 - Have a container image defined for every step.
 - Be aware of the different container bases you are using, and (for small or quick jobs) try to share them.
 - Don't be afraid to increase the number of jobs to a high number, as the flux scheduler on the MiniCluser is going to manage things.
 - Define your Kubernetes cluster size by exporting ``SNAKEMAKE_FLUX_OPERATOR_MAX_NODES``

This is a fairly simple workflow, meaning that it uses a single container and is able to run
all the jobs on the same MiniCluster.  