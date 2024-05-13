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

The Snakefile syntax obeys the following grammar, given in `extended Backus-Naur form (EBNF) <https://en.wikipedia.org/wiki/Extended_Backus%E2%80%93Naur_form>`_.

.. code-block:: text

    snakemake    = statement | rule | include | workdir | module | configfile | container
    rule         = "rule" (identifier | "") ":" ruleparams
    include      = "include:" stringliteral
    workdir      = "workdir:" stringliteral
    module       = "module" identifier ":" moduleparams
    configfile   = "configfile" ":" stringliteral
    userule      = "use" "rule" (identifier | "*") "from" identifier ["as" identifier] ["with" ":" norunparams]
    ni           = NEWLINE INDENT
    norunparams  = [ni input] [ni output] [ni params] [ni message] [ni threads] [ni resources] [ni log] [ni conda] [ni container] [ni benchmark] [ni cache] [ni priority]
    ruleparams   = norunparams [ni (run | shell | script | notebook)] NEWLINE snakemake
    input        = "input" ":" parameter_list
    output       = "output" ":" parameter_list
    params       = "params" ":" parameter_list
    log          = "log" ":" parameter_list
    benchmark    = "benchmark" ":" statement
    cache        = "cache" ":" bool
    message      = "message" ":" stringliteral
    threads      = "threads" ":" integer
    priority     = "priority" ":" integer
    resources    = "resources" ":" parameter_list
    version      = "version" ":" statement
    conda        = "conda" ":" stringliteral
    container    = "container" ":" stringliteral
    run          = "run" ":" ni statement
    shell        = "shell" ":" stringliteral
    script       = "script" ":" stringliteral
    notebook     = "notebook" ":" stringliteral
    moduleparams = [ni snakefile] [ni metawrapper] [ni config] [ni skipval]
    snakefile    = "snakefile" ":" stringliteral
    metawrapper  = "meta_wrapper" ":" stringliteral
    config       = "config" ":" stringliteral
    skipval      = "skip_validation" ":" stringliteral
    

while all not defined non-terminals map to their Python equivalents.

.. _snakefiles-depend_version:

Depend on a Minimum Snakemake Version
-------------------------------------

From Snakemake 3.2 on, if your workflow depends on a minimum Snakemake version, you can easily ensure that at least this version is installed via

.. code-block:: python

    from snakemake.utils import min_version

    min_version("3.2")

given that your minimum required version of Snakemake is 3.2. The statement will raise a WorkflowError (and therefore abort the workflow execution) if the version is not met.
