
.. _tutorial-flux:

Flux Tutorial
-------------

.. _Snakemake: http://snakemake.readthedocs.io
.. _Snakemake Remotes: https://snakemake.readthedocs.io/en/stable/snakefiles/remote_files.html
.. _Python: https://www.python.org/


Setup
:::::

To go through this tutorial, you need the following software installed:

- Docker


`Flux-framework <https://flux-framework.org/>`_ is a flexible resource scheduler that can work on both high performance computing systems and cloud (e.g., Kubernetes).
Since it is more modern (e.g., has an official Python API) we define it under a cloud resource. For this example, we will show you how to set up a "single node" local flux container to interact with snakemake. You can use the `Dockerfile in examples/flux <https://github.com/snakemake/snakemake/blob/main/examples/flux/Dockerfile>`_ that will provide a container with flux and snakemake        
Note that we install from source and bind to ``/home/fluxuser/snakemake`` with the intention of being able to develop (if desired).
First, build the container:

.. code-block:: console

    $ docker build -t flux-snake .

And then you can run the container with or without any such bind:

.. code-block:: console

    $ docker run -it --rm flux-snake 

Once you shelled into the container, you can view and start a flux instance:

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

See the `flux documentation <https://flux-framework.readthedocs.io/en/latest/quickstart.html#docker-recommended-for-quick-single-node-deployments>`_
for more detail. For now, let's try interacting with flux via snakemake via the `Flux Python Bindings <https://flux-framework.readthedocs.io/projects/flux-workflow-examples/en/latest/job-submit-api/README.html>`_.

The code for this example is provided in  (`examples/flux <https://github.com/snakemake/snakemake/tree/main/examples/flux>`_)
