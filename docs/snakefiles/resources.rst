.. _snakefiles-resources:

=========
Resources
=========

In addition to threads, a rule can use arbitrary user-defined resources by specifying them with the resources-keyword:

.. code-block:: python

    rule:
        input:     ...
        output:    ...
        resources: gpu=1
        shell: "..."

If limits for the resources are given via the command line, e.g.

.. code-block:: console

    $ snakemake --resources gpu=2

the scheduler will ensure that the given resources are not exceeded by running jobs.
If no limits are given, the resources are ignored.
Apart from making Snakemake aware of hybrid-computing architectures (e.g. with a limited number of additional devices like GPUs) this allows to control scheduling in various ways, e.g. to limit IO-heavy jobs by assigning an artificial IO-resource to them and limiting it via the ``--resources`` flag.
Resources must be ``int`` values.

Starting from version 3.7, resources can also be callables that return ``int`` values.
