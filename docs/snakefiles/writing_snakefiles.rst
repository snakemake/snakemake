.. _user_manual-writing_snakefiles:

=================
Writing Workflows
=================

In Snakemake, workflows are specified as Snakefiles.
Inspired by GNU Make, a Snakefile contains rules that denote how to create output files from input files.
Dependencies between rules are handled implicitly, by matching filenames of input files against output files.
Thereby wildcards can be used to write general rules.

.. _snakefiles-grammar:

-------
Grammar
-------

The Snakefile syntax obeys the following grammar, given in extended Backus-Naur form (EBNF)

.. code-block:: text

    snakemake  = statement | rule | include | workdir
    rule       = "rule" (identifier | "") ":" ruleparams
    include    = "include:" stringliteral
    workdir    = "workdir:" stringliteral
    ni         = NEWLINE INDENT
    ruleparams = [ni input] [ni output] [ni params] [ni message] [ni threads] [ni (run | shell)] NEWLINE snakemake
    input      = "input" ":" parameter_list
    output     = "output" ":" parameter_list
    params     = "params" ":" parameter_list
    log        = "log" ":" parameter_list
    benchmark  = "benchmark" ":" statement
    cache      = "cache" ":" bool
    message    = "message" ":" stringliteral
    threads    = "threads" ":" integer
    resources  = "resources" ":" parameter_list
    version    = "version" ":" statement
    run        = "run" ":" ni statement
    shell      = "shell" ":" stringliteral

while all not defined non-terminals map to their Python equivalents.

.. _snakefiles-depend_version:

Depend on a Minimum Snakemake Version
-------------------------------------

From Snakemake 3.2 on, if your workflow depends on a minimum Snakemake version, you can easily ensure that at least this version is installed via

.. code-block:: python

    from snakemake.utils import min_version

    min_version("3.2")

given that your minimum required version of Snakemake is 3.2. The statement will raise a WorkflowError (and therefore abort the workflow execution) if the version is not met.


.. _snakefiles-best_practices:

--------------
Best practices
--------------

Snakemake (>=5.11) comes with a code quality checker (a so called linter), that analyzes your workflow and highlights issues that should be solved in order to follow best practices, achieve maximum readability, and reproducibility.
The linter can be invoked with 

.. code-block:: bash

    snakemake --lint

given that a ``Snakefile`` or ``workflow/Snakefile`` is accessible from your working directory.
It is **highly recommended** to run the linter before publishing any workflow, asking questions on Stack Overflow or filing issues on Github.