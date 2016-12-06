.. _snakefiles-wrappers:

========
Wrappers
========

With Snakemake 3.5.5, the wrapper directive is introduced (experimental).
This directive allows to have re-usable wrapper scripts around e.g. command line tools. In contrast to modularization strategies like ``include`` or subworkflows, the wrapper directive allows to re-wire the DAG of jobs.
For example

.. code-block:: python

    rule samtools_sort:
        input:
            "mapped/{sample}.bam"
        output:
            "mapped/{sample}.sorted.bam"
        params:
            "-m 4G"
        threads: 8
        wrapper:
            "0.0.8/bio/samtools_sort"

Refers to the wrapper ``"0.0.8/bio/samtools_sort"`` to create the output from the input.
Snakemake will automatically download the wrapper from the [Snakemake Wrapper Repository](https://bitbucket.org/snakemake/snakemake-wrappers).
Thereby, 0.0.8 can be replaced with the git version tag you want to use, or a commit id (see [here](https://bitbucket.org/snakemake/snakemake-wrappers/commits/)).
This ensures reproducibility since changes in the wrapper implementation won't be propagated automatically to your workflow.
Alternatively, e.g., for development, the wrapper directive can also point to full URLs, including URLs to local files with ``file://``.
Examples for each wrapper can be found in the READMEs located in the wrapper subdirectories at the `Snakemake Wrapper Repository <https://bitbucket.org/snakemake/snakemake-wrappers>`_.

The `Snakemake Wrapper Repository <https://bitbucket.org/snakemake/snakemake-wrappers>`_ is meant as a collaborative project and pull requests are very welcome.
