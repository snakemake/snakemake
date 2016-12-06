.. _snakefiles-grammar:

=======
Grammar
=======

The snakefile syntax obeys the following grammar, given in extended Backus-Naur form (EBNF)

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
    message    = "message" ":" stringliteral
    threads    = "threads" ":" integer
    resources  = "resources" ":" parameter_list
    version    = "version" ":" statement
    run        = "run" ":" ni statement
    shell      = "shell" ":" stringliteral

while all not defined non-terminals map to their Python equivalents.