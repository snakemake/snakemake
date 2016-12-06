.. _snakefiles_cluster_configuration:

=====================
Cluster Configuration
=====================

Since version 3.3., Snakemake supports a separate configuration file for execution on a cluster. A cluster config file allows you to specify cluster submission parameters outside the Snakefile. The cluster config is a JSON- or YAML-formatted file that contains objects that match names of rules in the Snakefile. The parameters in the cluster config are then accessed by the cluster.* wildcard when you are submitting jobs. For example, say that you have the following Snakefile:

.. code-block: python

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

Any string in the cluster configuration can be formatted in the same way as shell commands, e.g. `"{rule}.{wildcards.sample}"` is formatted to `"a.xy"` if the rulename is `a` and the wildcard value is `xy`.
Here "__default__" is a special object that specifies default parameters, these will be inherited by the other configuration objects. The "compute1" object here changes the "time" parameter, but keeps the other parameters from "default". The rule "compute2" does not have any configuration, and will therefore use the default configuration. You can then run the Snakefile with the following command on a SLURM system.

.. code-block:: console

    $ snakemake -j 999 --cluster-config cluster.json --cluster "sbatch -A {cluster.account} -p {cluster.partition} -n {cluster.n}  -t {cluster.time}"
