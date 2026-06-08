.. _snakefiles-rules:

====================
Snakefiles and Rules
====================


.. toctree::
   :maxdepth: 1
   :hidden:

   rules/wildcards
   rules/executing_code
   rules/input_output
   rules/dependencies
   rules/reuse_composition
   rules/execution
   rules/directives
   rules/helpers


In Snakemake, workflows are specified as Snakefiles.
Inspired by GNU Make, a Snakefile contains rules that denote how to create output files from input files.
Dependencies between rules are handled implicitly, by matching filenames of input files against output files.
Thereby wildcards can be used to write general rules.


A Snakemake workflow defines a data analysis in terms of rules that are specified in the Snakefile.
Most commonly, rules consist of a name, input files, output files, and a shell command to generate the output from the input:

.. code-block:: python

    rule myrule:
        input:
            "path/to/inputfile",
            "path/to/other/inputfile",
        output:
            "path/to/outputfile",
            "path/to/another/outputfile",
        shell:
            "somecommand {input} {output}"

However, rules can be much more complex, may use :ref:`plain python <snakefiles-plain-python-rules>` or :ref:`boilerplate-free scripting in various languages <snakefiles-external_scripts>`, can contain :ref:`snakefiles-wildcards`, define :ref:`non-file parameters <snakefiles-params>`, :ref:`log files <snakefiles-log>` and many more, see below.

Inside the shell command, all local and global variables, especially input and output files can be accessed via their names in the `python format minilanguage <https://docs.python.org/py3k/library/string.html#formatspec>`_.
Here, input and output (and in general any list or tuple) automatically evaluate to a space-separated list of files (i.e. ``path/to/inputfile path/to/other/inputfile``).
From Snakemake 3.8.0 on, adding the special formatting instruction ``:q`` (e.g. ``"somecommand {input:q} {output:q}"``) will let Snakemake quote each of the list or tuple elements that contains whitespace.

.. note::

    Note that any placeholders in the shell command (like ``{input}``) are always evaluated and replaced
    when the corresponding job is executed, even if they are occurring inside a comment.
    To avoid evaluation and replacement, you have to mask the braces by doubling them,
    i.e. ``{{input}}``.

By default shell commands will be invoked with ``bash`` shell in the so-called  `strict mode <http://redsymbol.net/articles/unofficial-bash-strict-mode/>`_ (unless the workflow specifies something else, see :ref:`shell_settings`).


.. _snakefiles-depend_version:

Depend on a Minimum Snakemake Version
-------------------------------------

From Snakemake 3.2 on, if your workflow depends on a minimum Snakemake version, you can easily ensure that at least this version is installed via

.. code-block:: python

    from snakemake.utils import min_version

    min_version("3.2")

given that your minimum required version of Snakemake is 3.2. The statement will raise a WorkflowError (and therefore abort the workflow execution) if the version is not met.
