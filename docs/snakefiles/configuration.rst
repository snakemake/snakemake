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
