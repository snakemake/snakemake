.. _snakefiles-includes:

========
Includes
========

Another snakefile with all its rules can be included into the current:

.. code-block:

    include: "path/to/other/snakefile"

The default target rule (often called the "all"-rule), won't be affected by the include. I.e. it will always be the first rule in your Snakefile, no matter how many includes you have above your first rule.
From version 3.2 on, includes are relative to the directory of the Snakefile in which they occur. For example, if above Snakefile resides in the directory ``my/dir``, then Snakemake will search for the include at ``my/dir/path/to/other/snakefile``, regardless of the working directory.
