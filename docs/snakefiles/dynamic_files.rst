.. _snakefiles-dynamic_files:

=============
Dynamic Files
=============

Since version 1.2.3, Snakemake provides experimental support for dynamic files.
Dynamic files can be used whenever one has a rule, for which the number of output files is unknown before the rule was executed.
This is useful for example with cetain clustering algorithms:

.. code-block:: python

    rule cluster:
        input: "afile.csv"
        output: dynamic("{clusterid}.cluster.csv")
        run: ...

Now the results of the rule can be used in Snakemake although it does not know how many files will be present before executing the rule `cluster`, e.g. by:


.. code-block:: python

    rule all:
        input: dynamic("{clusterid}.cluster.plot.pdf")

    rule plot:
        input: "{clusterid}.cluster.csv"
        output: "{clusterid}.cluster.plot.pdf"
        run: ...

Here, Snakemake determines the input files for the rule `all` after the rule `cluster` was executed, and then dynamically inserts jobs of the rule `plot` into the DAG to create the desired plots.

Note that dynamic file support is still experimental. Especially, using more than one wildcard within dynamic files can introduce various problems.
Before using dynamic files, think about alternative, static solutions, where you know beforehand how many output files your rule will produce. In four years and hundreds of workflows, I needed dynamic files only once.
