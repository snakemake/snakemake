================
Helper functions
================

.. _snakefiles-input_helpers:

Input Helpers
--------------------------

Snakemake provides a number of helpers that can be used to define rules and drastically simplify over using
:ref:`input functions <snakefiles-input_functions>` or :ref:`plain python expressions <snakefiles_aggregation>`.
Below, we will first start with describing two basic helper functions for specifying aggregations and multiple output files.
Afterwards, we will further show a set of semantic helper functions should increase readability and simplify code (see :ref:`snakefiles-semantic-helpers`).

.. _snakefiles_expand:

The expand function
~~~~~~~~~~~~~~~~~~~

Instead of specifying input files via a Python list comprehension, Snakemake offers a helper function ``expand()``.

.. code-block:: python

    rule aggregate:
        input:
            expand("{dataset}/a.txt", dataset=DATASETS)
        output:
            "aggregated.txt"
        shell:
            ...


Note that *dataset* is NOT a wildcard here because it is resolved by Snakemake due to the ``expand`` statement.
The ``expand`` function also allows us to combine different variables, e.g.

.. code-block:: python

    rule aggregate:
        input:
            expand("{dataset}/a.{ext}", dataset=DATASETS, ext=FORMATS)
        output:
            "aggregated.txt"
        shell:
            ...

If ``FORMATS=["txt", "csv"]`` contains a list of desired output formats then expand will automatically combine any dataset with any of these extensions.
Furthermore, the first argument can also be a list of strings. In that case, the transformation is applied to all elements of the list. E.g.

.. code-block:: python

    expand(["{dataset}/a.{ext}", "{dataset}/b.{ext}"], dataset=DATASETS, ext=FORMATS)

leads to

.. code-block:: python

    ["ds1/a.txt", "ds1/b.txt", "ds2/a.txt", "ds2/b.txt", "ds1/a.csv", "ds1/b.csv", "ds2/a.csv", "ds2/b.csv"]

Per default, ``expand`` uses the python itertools function ``product`` that yields all combinations of the provided wildcard values. However by inserting a second positional argument this can be replaced by any combinatoric function, e.g. ``zip``:

.. code-block:: python

    expand(["{dataset}/a.{ext}", "{dataset}/b.{ext}"], zip, dataset=DATASETS, ext=FORMATS)

leads to

.. code-block:: python

    ["ds1/a.txt", "ds1/b.txt", "ds2/a.csv", "ds2/b.csv"]

You can also mask a wildcard expression in ``expand`` such that it will be kept, e.g.

.. code-block:: python

    expand("{{dataset}}/a.{ext}", ext=FORMATS)

will create strings with all values for ext but starting with the wildcard ``"{dataset}"``.

Finally, argument values passed to ``expand`` can also be functions or lists of functions if the return value of ``expand`` or ``expand`` itself is used within ``input``, or ``params``.
Depending on the context, that function has to accept the same arguments as functions for ``input`` (see :ref:`snakefiles-input_functions`) or functions for ``params`` (see :ref:`snakefiles-params`).
If that is the case, ``expand`` returns a function again, the evaluation of which is deferred to the point in time when the wildcards of the respective job are known.


.. _snakefiles-multiext:

The multiext function
~~~~~~~~~~~~~~~~~~~~~

``multiext`` provides a simplified variant of ``expand`` that allows us to define a set of output or input files that just differ by their extension:


.. code-block:: python

    rule plot:
        input:
            ...
        output:
            multiext("some/plot", ".pdf", ".svg", ".png")
        shell:
            ...

The effect is the same as if you would write ``expand("some/plot{ext}", ext=[".pdf", ".svg", ".png"])``, however, using a simpler syntax.
Moreover, defining output with ``multiext`` is the only way to use :ref:`between workflow caching <caching>` for rules with multiple output files.
It's also possible to get named input/output files in the following way:

.. code-block:: python

    rule plot:
        input:
            ...
        output:
            multiext("some/plot", out1=".pdf", out2=".svg")
            "some_other_output"
            named_output="another_output"
        shell:
            """
            somecommand > {output.out1}
            othercommand > {output.out2}
            anothercommand > {output[2]}
            finalcommand > {output.named_output}
            """

Do note that all the multiext extensions should be named, or all of them should be unnamed (not both).
Additionally, if additional input/output statements are given, multiext should be treated as positional arguments (before other named input/output files).



.. _snakefiles-semantic-helpers:

Semantic helpers
~~~~~~~~~~~~~~~~

The collect function
""""""""""""""""""""

The ``collect`` function is an alias for the ``expand`` function with exactly the same behavior.
It can be used to express more explicitly that a rule collects a set of files from upstream jobs.

The lookup function
"""""""""""""""""""

The ``lookup`` function can be used to look up a value in a python mapping (e.g. a ``dict``) or a `pandas dataframe or series <https://pandas.pydata.org>`_.
It is especially useful for looking up information based on wildcard values.
The ``lookup`` function has the signature

.. code-block:: python

    lookup(
        dpath: Optional[str | Callable] = None,
        query: Optional[str | Callable] = None,
        cols: Optional[List[str]] = None,
        is_nrows: Optional[int],
        within=None,
        default=NODEFAULT
    )

The required ``within`` parameter takes either a python mapping, a pandas dataframe, or a pandas series.
For the former case, it expects the ``dpath`` argument, for the latter two cases, it expects the ``query`` argument to be given.

In case of a pandas dataframe,
the query parameter is passed to `DataFrame.query() <https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.query.html>`_.
If the query results in multiple rows, the result is returned as a list of
named tuples with the column names as attributes.
If the query results in a single row, the result is returned as a single
named tuple with the column names as attributes.
If the query or dpath parameter is given a function, the function will be evaluated with wildcards passed as the first argument.
If the dpath is not found or the query returns no matching rows, the ``default`` fallback
value is returned if provided. Otherwise, a ``LookupError`` is raised (for dpath) or an empty list is returned (for query).
Note: ``None`` is also a valid default value.

In both cases (``dpath`` and ``query``), the result can be used by the ``expand`` or ``collect`` function,
e.g.

.. code-block:: python

    collect("results/{item.sample}.txt", item=lookup(query="someval > 2", within=samples))

Here, we take the file ``"results/{item.sample}.txt"`` with ``{item.sample}`` being replaced by the
sample names that occur in all rows of the dataframe ``samples`` where the value of the ``someval`` column is greater than 2.

Since the result, in any case, also evaluates to True if it is not empty
when interpreted as a boolean by Python, it can also be used as a condition
for the :ref:`branch function <snakefiles-branch-function>`, e.g.

.. code-block:: python

    branch(lookup(query="sample == '{sample}' & someval > 2", within=samples), then="foo", otherwise="bar")

In case your dataframe has an index, you can also access the index within the
query, e.g. for faster, constant time lookups:

.. code-block:: python

    lookup(query="index.loc[{sample}]", within=samples)

Further, it is possible to constrain the output to a list of columns, e.g.

.. code-block:: python

    lookup(query="sample == '{sample}'", within=samples, cols=["somecolumn"])

or to a single column, e.g.

.. code-block:: python

    lookup(query="sample == '{sample}'", within=samples, cols="somecolumn")

In the latter case, just a list of items in that column is returned (e.g. ``["a", "b", "c"]``).

The argument ``is_nrows`` allows to test for a given number of rows in the queried dataframe.
If it is used, lookup just returns a boolean value indicating whether the number of rows in the queried dataframe matches the given number:

.. code-block:: python

    lookup(query="sample == '{sample}'", within=samples, is_nrows=5)

In case of a **pandas series**, the series is converted into a dataframe via
Series.to_frame() and the same logic as for a dataframe is applied.

In case of a **python mapping**, the dpath parameter is passed to dpath.values()
(see https://github.com/dpath-maintainers/dpath-python).

``query``, ``dpath``, and ``cols`` may contain wildcards (e.g. ``{sample}``).
In that case, this function returns an :ref:`input function <snakefiles-input_functions>` which takes
wildcards as its only argument and will be evaluated by Snakemake
once the wildcard values are known if the lookup is used within an input file statement.

In addition to wildcard values, dpath, query and cols may refer via the same syntax
to auxiliary namespace arguments given to the lookup function, e.g.

.. code-block:: python

    lookup(
        query="cell_type == '{sample.cell_type}'",
        within=samples,
        sample=lookup("sample == '{sample}'", within=samples)
    )

This way, one can e.g. pass additional variables or chain lookups into more complex queries.

.. _snakefiles-branch-function:

The branch function
"""""""""""""""""""

The ``branch`` function allows to choose different input files based on a given conditional.
It has the signature

.. code-block:: python

    branch(
        condition: Union[Callable, bool],
        then: Optional[Union[str, list[str], Callable]] = None,
        otherwise: Optional[Union[str, list[str], Callable]] = None,
        cases: Optional[Mapping] = None
    )

The ``condition`` argument has to be either a function or an expression that can be evaluated as a ``bool`` (which is virtually everything in Python).
If it is a function, it has to take wildcards as its only parameter.
Similarly, ``then``, ``otherwise`` and the values of the ``cases`` mapping (e.g. a python ``dict``) can be such functions.

If any such function is given to any of those arguments, this function returns a derived
input function that will be evaluated once the wildcards are known (e.g. when used in the context of an input definition) (see :ref:`snakefiles-input_functions`).

If ``then`` and optionally ``otherwise`` are specified, it does the following:
If the ``condition`` is (or evaluates to) ``True``, return the value
of the ``then`` parameter. Otherwise, return the value of the ``otherwise`` parameter.

If ``cases`` is specified, it does the following:
Retrieve the value of the cases mapping using the return value of the condition
(if it is a function), or the condition value itself as a key.

An example of using ``branch`` in combination with ``lookup`` from a ``config`` dictionary can look as follows:

.. code-block:: python

    branch(
        lookup(dpath="tools/sometool", within=config),
        then="results/sometool/{dataset}.txt",
        otherwise="results/someresult/{dataset}.txt"
    )

Here, the semantic is as follows:
If the lookup returns ``True``, the input is ``results/sometool/{dataset}.txt``, otherwise it is ``results/someresult/{dataset}.txt``.

Given that ``condition`` can be a function, if this is used in the context of a rule definition and the usage of the tool ``sometool`` depends on some wildcard values,
one can also pass a function name instead of a boolean value to the branch function (using it as an input function).

.. code-block:: python

    def use_sometool(wildcards):
        # determine whether the tool shall be used based on the wildcard values.
        ...

    rule a:
        input:
            branch(
                use_sometool,
                then="results/sometool/{dataset}.txt",
                otherwise="results/someresult/{dataset}.txt"
            )

Above, the semantic is as follows:
If ``use_sometool`` returns ``True`` for the given wildcard values, the input is ``results/sometool/{dataset}.txt``, otherwise it is ``results/someresult/{dataset}.txt``.

An example for using the cases argument could look as follows:

.. code-block:: python

    branch(
        lookup(dpath="tool/to/use", within=config),
        cases={
            "sometool": "results/sometool/{dataset}.txt",
            "someothertool": "results/someothertool/{dataset}.txt"
        }
    )

The evaluate function
"""""""""""""""""""""

The ``evaluate`` function allows to quickly evaluate a Python expression that contains wildcard values.
It has the signature ``evaluate(expr: str)``.
Within the expression one can specify wildcards via the usual syntax, e.g. ``{sample}``.
Upon evaluation, the wildcards are replaced by their values as strings and the expression is evaluated as Python code with access to any global variables defined in the workflow.
Consider the following example:

.. code-block:: python

    rule a:
    input:
        branch(evaluate("{sample} == '100'"), then="a/{sample}.txt", otherwise="b/{sample}.txt"),
    output:
        "c/{sample}.txt",
    shell:
        ...

The semantic is as follows:
If the sample wildcard is ``100``, the input is ``a/100.txt``, otherwise it is ``b/100.txt``.

.. _snakefiles-semantic-helpers-exists:

The exists function
"""""""""""""""""""

The ``exists`` function allows to check whether a file exists, while properly considering remote storage settings provided to Snakemake.
For example, if Snakemake has been configured to consider all input and output files to be located in an S3 bucket, ``exists`` will check whether the file exists in the S3 bucket.
It has the signature ``exists(path)``, with ``path`` being the path to a file or directory, or an explicit :ref:`storage object <storage-support>`.
The function returns ``True`` if the file exists, and ``False`` otherwise.
It can for example be used to condition some behavior in the workflow on the existence of a file **before** the workflow is executed:

.. code-block:: python

    rule all:
        input:
            # only expect the output if test.txt is present before workflow execution
            "out.txt" if exists("test.txt") else [],

    rule b:
        input:
            "test.txt"
        output:
            "out.txt"
        shell:
            "cp {input} {output}"

.. _snakefiles-semantic-helpers-parse-input:

The parse_input function
""""""""""""""""""""""""

The ``parse_input`` function allows to parse an input file and return a value.
It has the signature ``parse_input(input_item, parser, kwargs)``, with ``input_item`` being the key of an input file, ``parser`` being a callable to extract the desired information, and ``kwargs`` extra arguments passed to the parser.
The function will return the extracted value and this can, for example, be used as a parameter.

.. code-block:: python

    rule a:
	input:
	    samples="samples.tsv",
        output:
            "samples.id",
        params:
            id=parse_input(input.samples, parser=extract_id)
        shell:
            "echo {params.id} > {output}"


.. _snakefiles-semantic-helpers-extract-checksum:

The extract_checksum function
"""""""""""""""""""""""""""""

The ``extract_checksum`` function parses an input file and returns the checksum of the given file.
It has the signature ``extract_checksum(infile, file)``, with ``infile`` being the input file, and ``file`` the filename to search for.
The function will return the checksum of ``file`` present in ``infile``.

.. code-block:: python

    rule a:
	input:
	    checksum="samples.md5",
        output:
            tsv="{a}.tsv",
        params:
            checksum=parse_input(input.checksum, parser=extract_checksum, file=output.tsv)
        shell:
            "echo {params.checksum} > {output}"

.. _snakefiles-rule-item-access:

Rule item access helpers
""""""""""""""""""""""""

Via functions (e.g. for :ref:`snakefiles-params` or :ref:`snakefiles-resources`) it is possible to access other items of the same rule in a deferred way, at the point in time when they are actually known.
For this, functions like

.. code-block:: python

    def get_file_foo_from_input(wildcards, input):
        return input.foo

can be written.
If such a function is passed to e.g. a params or resource statement, Snakemake knows that this resource shall be evaluated by passing the input files in addition to the wildcards (which are always required as first argument for any such function).
To simplify such logic for certain situations, Snakemake provides globally available objects
``input``, ``output``, ``resources``, and ``threads`` that can be used to replace the corresponding function definitions.
For example, the global ``input.foo`` (not the one inside above function, which returns its value from the ``input`` argument of the function, which in turn is a concrete file path) returns a function that is equivalent to ``get_file_foo_from_input`` (the function above).
Using these objects makes most sense inside of a rule definition.
For example, it can be used to access a subpath of an input or output file or directory, see :ref:`snakefiles-subpath`.
For example, we could write

.. code-block:: python

    rule a:
        input:
            foo="results/something/foo.txt"
        output:
            "results/something-else/out.txt"
        params:
            directory=subpath(input.foo, parent=True)
        shell:
            "somecommand {params.directory} {output}"

.. _snakefiles-subpath:

Sub-path access
"""""""""""""""

In some cases, it is useful to access a sub-path of an input or output file or directory.
For this purpose, Snakemake provides the ``subpath`` function.
It has the signature ``subpath(path_or_func, strip_suffix=None, basename=False, parent=False, ancestor=None)``.
If a path is given as first argument (of type ``str`` or ``pathlib.Path``), the function directly returns the sub-path of the given path.
Thereby, the sub-path is determined depending on the other arguments.

If a ``str`` is given to ``strip_suffix``, this suffix is stripped from the path before determining the sub-path (a ``ValueError`` error is thrown if the path does not have the suffix).

.. code-block:: python

    subpath("results/test.txt", strip_suffix=".txt") # returns "results/test"

If ``basename`` is set to ``True``, the basename of the path is returned (e.g. ``test.txt`` in case the path is ``results/test.txt``).

.. code-block:: python

    subpath("results/test.txt", basename=True) # returns "test.txt"

If ``parent`` is set to ``True``, the parent directory of the path is returned (e.g. ``results`` in case the path is ``results/test.txt``).

.. code-block:: python

    subpath("results/test.txt", parent=True) # returns "results"

If ``ancestor`` is set to an integer greater than 0, the ancestor directory at the given level is returned (e.g. ``results`` in case the path is ``results/foo/test.txt`` and ``ancestor=2``).

.. code-block:: python

    subpath("results/foo/test.txt", ancestor=2) # returns "results"

The arguments ``basename``, ``parent``, and ``ancestor`` are mutually exclusive.

The ``subpath`` function can be very handy in combination with :ref:`Snakemake's rule item access helpers <snakefiles-rule-item-access>`, e.g.

.. code-block:: python

    rule a:
        input:
            "results/something/foo.txt"
        output:
            foo="results/something-else/out.txt"
        params:
            basename=subpath(output.foo, basename=True),
            outdir=subpath(output.foo, parent=True)
        shell:
            "somecommand {input} --name {params.basename} --outdir {params.outdir}"


.. _snakefiles-flatten:

flatten
"""""""
When selecting input files, sometimes you might end up with an irregular list of lists. To flatten in, you can use:

.. code-block:: python

    flatten([1, "a", [2,"b"], ["c","d",["e", 3]]]) # returns ["1", "a", "2", "b", "c", "d", "e", "3"]



    .. _snakefiles-pathvars:

    Path variables
    ~~~~~~~~~~~~~~

    Certain components in input and output file paths tend to reoccur across many rules.
    Via so-called **pathvars**, Snakemake allows to define such components globally, make them configurable via the config file, and change them per module or even per rule.
    Apart from saving boilerplate code, pathvars can be used to make modules intended for reuse in multiple contexts more flexible.
    Pathvars can be used as generic placeholders for their actual values inside of input, output, log, and benchmark paths, using angle brackets, e.g. ``<results>``.
    They behave similarly to `Python string interpolation <https://docs.python.org/3/tutorial/inputoutput.html#formatted-string-literals>`__ but only allow predefined placeholders, with precedence and configuration handled by Snakemake.

    Pathvar usage
    """""""""""""

    An example rule using pathvars is the following:

    .. code-block:: python

        rule somerule:
            input:
                "<results>/something/{sample}.txt"
            output:
                "<results>/processed/{sample}.txt"
            shell:
                "somecommand {input} {output}"

    Pathvars are resolved when the rule is parsed (before wildcard resolution and DAG construction).
    The values of pathvars can thereby even contain wildcards themselves.

    Pathvar defaults
    """"""""""""""""

    By default, Snakemake offers the pathvars ``results``, ``resources``, ``logs``, ``benchmarks``.
    Each of them is set to its respective name (i.e. the output file ``"<results>/processed/{sample}.txt"`` will be interpreted as ``"results/processed/{sample}.txt"``).

    Pathvar definition
    """"""""""""""""""

    Beyond the defaults, it is possible to define additional pathvars or customize the default definitions.
    This can happen in multiple ways, with the following precedence (from highest to lowest):

    1. For individual rules, via the ``pathvars`` keyword.
    2. For :ref:`module <snakefiles-modules>` config via the ``pathvars`` key in a config dict explicitly passed to the module (this applies recursively to nested modules).
    3. For :ref:`modules <snakefiles-modules>`, via the ``pathvars`` keyword to the ``module`` directive (this applies recursively to nested modules).
    4. Globally, via the ``pathvars`` key in the config file or passed to the ``--config`` command line arguments.
    5. Globally, via the ``pathvars`` keyword at the top level of the Snakefile.

    Thereby, if two definitions share the same precedence, the last one wins.

    Apart from :ref:`module pathvars <snakefiles-modules-pathvars>`, the most common way is to define them globally via the ``pathvars`` keyword:

    .. code-block:: python

        pathvars:
            per="{sample}"

    Above, we define a pathvar ``per`` and set it to the value ``{sample}``, thus defining a wildcard.
    Such a pattern can be helpful if you write a workflow that shall be reused as a module in various different ways within other workflows (e.g. thereby processing different items, samples, or something else).

    In order to overwrite pathvars for individual rules, they can be specified via the ``pathvars`` keyword inside a rule:

    .. code-block:: python

        rule somerule:
            input:
                "<results>/something/<per>.txt"
            output:
                "<results>/processed/<per>.txt"
            pathvars:
                results="custom-folder"
            shell:
                "somecommand {input} {output}"

    This way, the pathvar ``<results>`` in the input and output path would be replaced with ``custom-folder`` just for this rule.
    Per rule pathvar definition can also happen in combination with :ref:`rule inheritance <snakefiles-rule-inheritance>`.
    This allows to quickly write and reuse rules with generic input or output files.

    Finally, it is possible overwrite pathvars via the workflow configuration (configfile or ``--config``).
    For this purpose, it is possible to define a key ``pathvars`` in the config, with a mapping between pathvars and their values below, e.g.

    .. code-block:: yaml

        pathvars:
            results: example-folder

    Note that defining pathvars in the config should be considered a rare, discouraged, and advanced use case, since the user must know the workflow's internal pathvar expectations.
    Workflow authors can explicitly forbid the modification of particular pathvars via :ref:`config file schemas and validation <snakefiles_config_validation>`.
