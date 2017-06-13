.. _snakefiles_configuration:

=============
Configuration
=============

Snakemake allows you to use configuration files for making your workflows more flexible and also for abstracting away direct dependencies to a fixed HPC cluster scheduler.


.. _snakefiles_standard_configuration:

----------------------
Standard Configuration
----------------------

Snakemake directly supports the configuration of your workflow.
A configuration is provided as a JSON or YAML file and can be loaded with:

.. code-block:: python

    configfile: "path/to/config.json"

The config file can be used to define a dictionary of configuration parameters and their values.
In the workflow, the configuration is accessible via the global variable `config`, e.g.

.. code-block:: python

    rule all:
        input:
            expand("{sample}.{yourparam}.output.pdf", sample=config["samples"], param=config["yourparam"])

If the `configfile` statement is not used, the config variable provides an empty array.
In addition to the `configfile` statement, config values can be overwritten via the command line or the :ref:`api_reference_snakemake`, e.g.:

.. code-block:: console

    $ snakemake --config yourparam=1.5

Further, you can manually alter the config dictionary using any Python code **outside** of your rules. Changes made from within a rule won't be seen from other rules.
Finally, you can use the `--configfile` command line argument to overwrite values from the `configfile` statement.
Note that any values parsed into the `config` dictionary with any of above mechanisms are merged, i.e., all keys defined via a `configfile`
statement, or the `--configfile` and `--config` command line arguments will end up in the final `config` dictionary, but if two methods define the same key, command line
overwrites the `configfile` statement.

For adding config placeholders into a shell command, Python string formatting syntax requires you to leave out the quotes around the key name, like so:

.. code-block:: python

    shell:
        "mycommand {config[foo]} ..."


.. _snakefiles-cluster_configuration:

---------------------
Cluster Configuration
---------------------

Snakemake supports a separate configuration file for execution on a cluster.
A cluster config file allows you to specify cluster submission parameters outside the Snakefile.
The cluster config is a JSON- or YAML-formatted file that contains objects that match names of rules in the Snakefile.
The parameters in the cluster config are then accessed by the ``cluster.*`` wildcard when you are submitting jobs.
For example, say that you have the following Snakefile:

.. code-block:: python

    rule all:
        input: "input1.txt", "input2.txt"

    rule compute1:
        output: "input1.txt"
        shell: "touch input1.txt"

    rule compute2:
        output: "input2.txt"
        shell: "touch input2.txt"

This Snakefile can then be configured by a corresponding cluster config, say "cluster.json":


.. code-block:: json

    {
        "__default__" :
        {
            "account" : "my account",
            "time" : "00:15:00",
            "n" : 1,
            "partition" : "core"
        },
        "compute1" :
        {
            "time" : "00:20:00"
        }
    }

Any string in the cluster configuration can be formatted in the same way as shell commands, e.g. ``{rule}.{wildcards.sample}`` is formatted to ``a.xy`` if the rulename is ``a`` and the wildcard value is ``xy``.
Here ``__default__`` is a special object that specifies default parameters, these will be inherited by the other configuration objects. The ``compute1`` object here changes the ``time`` parameter, but keeps the other parameters from ``__default__``. The rule ``compute2`` does not have any configuration, and will therefore use the default configuration. You can then run the Snakefile with the following command on a SLURM system.

.. code-block:: console

    $ snakemake -j 999 --cluster-config cluster.json --cluster "sbatch -A {cluster.account} -p {cluster.partition} -n {cluster.n}  -t {cluster.time}"


For cluster systems using LSF/BSUB, a cluster config may look like this:

.. code-block:: json

    {
        "__default__" :
        {
            "queue"     : "medium_priority",
            "nCPUs"     : "16",
            "memory"    : 20000,
            "resources" : "\"select[mem>20000] rusage[mem=20000] span[hosts=1]\"",
            "name"      : "JOBNAME.{rule}.{wildcards}",
            "output"    : "logs/cluster/{rule}.{wildcards}.out",
            "error"     : "logs/cluster/{rule}.{wildcards}.err"
        },


        "trimming_PE" :
        {
            "memory"    : 30000,
            "resources" : "\"select[mem>30000] rusage[mem=30000] span[hosts=1]\"",
        }
    }

The advantage of this setup is that it is already pretty general by exploiting the wildcard possibilities that Snakemake provides via ``{rule}`` and ``{wildcards}``. So job names, output and error files all have reasonable and trackable default names, only the directies (``logs/cluster``) and job names (``JOBNAME``) have to adjusted accordingly.
If a rule named ``bamCoverage`` is executed with the wildcard ``basename = sample1``, for example, the output and error files will be ``bamCoverage.basename=sample1.out`` and ``bamCoverage.basename=sample1.err``, respectively.


---------------------------
Configure Working Directory
---------------------------

All paths in the snakefile are interpreted relative to the directory snakemake is executed in. This behaviour can be overridden by specifying a workdir in the snakefile:

.. code-block:: python

    workdir: "path/to/workdir"

Usually, it is preferred to only set the working directory via the command line, because above directive limits the portability of Snakemake workflows.
