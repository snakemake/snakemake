.. _snakefiles-grammar:

=======
Grammar
=======

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