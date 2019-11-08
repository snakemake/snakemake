==========================================================
Caching and reusing intermediate results between workflows
==========================================================

Within certain data analysis fields, there are certain intermediate results that reoccur in exactly the same way in many analysis.
For example, in bioinformatics, reference genomes or annotations are downloaded, and read mapping indexes are built.
Since such steps are independent of the actual data or measurements that are analyzed, but still computationally or timely expensive to conduct, it has been common practice to externalize their computation and assume the presence of the resulting files before execution of a workflow.

From version 5.8.0 on, Snakemake offers a way to keep those steps inside the actual analysis without requiring from redundant computations.
By hashing all steps, parameters, software stacks (in terms of conda environments or containers), and raw input required up to a certain step in a `blockchain <https://en.wikipedia.org/wiki/Blockchain>`_, Snakemake is able to recognize **before** the computation whether a certain result is already available in a central cache on the same system.

Such caching has to be explitly activated per rule, which can be done via the command line interface.
For example,

.. code-block:: console

    $ export SNAKEMAKE_OUTPUT_CACHE=/mnt/snakemake-cache/
    $ snakemake --cache download_reference create_index

would instruct Snakemake to cache and reuse the results of the rules ``download_reference``and ``create_index``.
The environment variable definition that happens in the first line (defining the location of the cache) should of course be done only once and system wide in reality.

Note that only rules with just a single output file are eligible for caching.
Also note that the rules need to retrieve all their parameters via the ``params`` directive (except input files).
It is not allowed to directly use ``wildcards``, ``config`` or any global variable in the shell command or script, because these are not captured in the hash (otherwise, reuse would be unnecessarily limited).

Also note that Snakemake will store everything in the cache as readable and writeable for all users on the system.
Hence, caching is not intended for private data, just for steps that deal with publicly available resources.

Finally, be aware that the implementation has to be considered **experimental** until this note is removed.