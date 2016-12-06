.. _snakefiles-input_functions:

===============
Input Functions
===============

Instead of specifying strings or lists of strings as input files, snakemake can also make use of functions that return single **or** lists of input files:

.. code-block:: python

    def myfunc(wildcards):
        return [... a list of input files depending on given wildcards ...]

    rule:
        input: myfunc
        output: "someoutput.{somewildcard}.txt"
        shell: "..."

The function has to accept a single argument that will be the wildcards object generated from the application of the rule to create some requested output files.
By this, rules can have entirely different input files (both in form and number) depending on the inferred wildcards. E.g. you can assign input files that appear in entirely different parts of your filesystem based on some wildcard value and a dictionary that maps the wildcard value to file paths.

Note that the function will be executed when the rule is evaluated and before the workflow actually starts to execute. Further note that using a function as input overrides the default mechanism of replacing wildcards with their values inferred from the output files. You have to take care of that yourself with the given wildcards object.

Finally, when implementing the input function, it is best practice to make sure that it can properly handle all possible wildcard values your rule can have.
In particular, input files should not be combined with very general rules that can be applied to create almost any file: Snakemake will try to apply the rule, and will report the exceptions of your input function as errors.
