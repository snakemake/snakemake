==========================
Rule Reuse and Composition
==========================

.. _snakefiles-includes:

--------
Includes
--------

Another Snakefile with all its rules can be included into the current:

.. code-block:: python

    include: "path/to/other/snakefile"

The default target rule (often called the ``all``-rule), won't be affected by the include.
I.e. it will always be the first rule in your Snakefile, no matter how many includes you have above your first rule.
Includes are relative to the directory of the Snakefile in which they occur.
For example, if above Snakefile resides in the directory ``my/dir``, then Snakemake will search for the include at ``my/dir/path/to/other/snakefile``, regardless of the working directory.


.. _snakefiles-modules:

-------
Modules
-------

With Snakemake 6.0 and later, it is possible to define external workflows as modules, from which rules can be used by explicitly "importing" them.

.. code-block:: python

    from snakemake.utils import min_version
    min_version("6.0")

    module other_workflow:
        snakefile:
            # here, plain paths, URLs and the special markers for code hosting providers (see below) are possible.
            "other_workflow/Snakefile"

    use rule * from other_workflow exclude ruleC as other_*

The ``module other_workflow:`` statement registers the external workflow as a module, by defining the path to the main snakefile of ``other_workflow``.
Here, plain paths, HTTP/HTTPS URLs and special markers for code hosting providers like Github or Gitlab are possible (see :ref:`snakefile-code-hosting-providers`).
The second statement, ``use rule * from other_workflow exclude ruleC as other_*``, declares all rules of that module to be used in the current one, except for ruleC.
Thereby, the ``as other_*`` at the end renames all those rules with a common prefix.
This can be handy to avoid rule name conflicts (note that rules from modules can otherwise overwrite rules from your current workflow or other modules).

.. note::

    The imported module cannot be named as `workflow`, which is a reserved name.

The module is evaluated in a separate namespace, and only the selected rules are added to the current workflow.
Non-rule Python statements inside the module are also evaluated in that separate namespace.
They are available in the module-defining workflow under the name of the module (e.g. here ``other_workflow.myfunction()`` would call the function ``myfunction`` that has been defined in the model, e.g. in ``other_workflow/Snakefile``).
Also note that this means that any Python variables and functions available in the module-defining namespace will **not** be visible from inside the module.
However, it is possible to pass information to the module using the ``config`` mechanism described in the following.

It is possible to overwrite the global config dictionary for the module, which is usually filled by the ``configfile`` statement (see :ref:`snakefiles_standard_configuration`):

.. code-block:: python

    from snakemake.utils import min_version
    min_version("6.0")

    configfile: "config/config.yaml"

    module other_workflow:
        # here, plain paths, URLs and the special markers for code hosting providers (see below) are possible.
        snakefile: "other_workflow/Snakefile"
        config: config["other-workflow"]

    use rule * from other_workflow as other_*

In this case, any ``configfile`` statements inside the module are ignored.
In addition, it is possible to skip any :ref:`validation <snakefiles_config_validation>` statements in the module, by specifying ``skip_validation: True`` in the module statement.
Moreover, one can automatically move all relative input and output files of a module into a dedicated folder by specifying ``prefix: "foo"`` in the module definition, e.g. any output file ``path/to/output.txt`` in the module would be stored under ``foo/path/to/output.txt`` instead.
This becomes particularly useful when combining multiple modules, see :ref:`use_with_modules`.
However, if you have some input files that come from outside the workflow, you can use the ``local`` flag so that their path is not modified (see :ref:`snakefiles-storage-local-files`).

Instead of using all rules, you can selectively import specific rules from modules.
These rules can also be modified during import via a final ``with:`` followed by a block that lists items to overwrite.
This behaves similarly to inheriting an existing rule within the current workflow, but with a ``from`` statement to declare the original module (see :ref:`snakefiles-rule-inheritance`).

.. code-block:: python

    from snakemake.utils import min_version
    min_version("6.0")

    module other_workflow:
        # here, plain paths, URLs and the special markers for code hosting providers (see below) are possible.
        snakefile: "other_workflow/Snakefile"
        config: config["other-workflow"]

    use rule * from other_workflow as other_*

    use rule some_task from other_workflow as other_some_task with:
        output:
            "results/some-result.txt"

When using the ``with:`` block, keyword arguments in ``params`` will be selectively replaced, while positional arguments are overwritten if provided.
All other properties (e.g., ``input``, ``output``, ``log``, ``params``, etc.) will be fully overwritten with the values specified in the block, except the actual execution step (``shell``, ``notebook``, ``script``, ``cwl``, or ``run``).

Note that the second use statement has to use the original rule name, not the one that has been prefixed with ``other_`` via the first use statement (there is no rule ``other_some_task`` in the module ``other_workflow``).

.. note::

    A rule cannot be overwritten under the same name, unless it was previously imported via `use rule * from ...` statement.
    This is the **only allowed scenario** where an existing rule name may be overwritten, and is provided for convenience when selectively customizing some rules without introducing new names.
    In such cases, the second statement uses the same final name as produced by the previous import (via the `as` clause).
    Importantly, once a rule has been modified in this way, it cannot be redefined or modified again under the same name, but you should import under different names to customize the same rule multiple times:

    .. code-block:: python

        use rule * from other_workflow as other_*

        use rule some_task from other_workflow as other_some_task with:
            output:
                "results/some-result.txt"

        use rule some_task from other_workflow as else_some_task with:
            output:
                "custom_output.txt"

    Once a rule has been modified this way under a given name, it **cannot** be redefined or modified again under the same name:

    .. code-block:: python

        use rule some_task from other_workflow as other_some_task with:
            output:
                "results/some-result.txt"

        use rule some_task from other_workflow as other_some_task with:
            threads: 1
        # Not allowed: "other_some_task" was already defined above.

    Similarly, if a `use rule * from ...` statement would result in a rule name that collides with a previously defined rule (regardless of its source), Snakemake will raise an error, and you should resolve the conflict by changing the import order or using a different `as` modifier:

    .. code-block:: python

        use rule some_task from other_workflow as else_some_task with:
            output:
                "custom_output.txt"

        use rule * from other_workflow as else_*
        # Will fail: "else_some_task" is already defined.

Of course, it is possible to combine the use of rules from multiple modules (see :ref:`use_with_modules`), and via modifying statements they can be rewired and reconfigured in an arbitrary way.

.. _snakefiles-rule-inheritance:

Rule inheritance
----------------

With Snakemake 6.0 and later, it is possible to inherit from previously defined rules, or in other words, reuse an existing rule in a modified way.
This works via the ``use rule`` statement that also allows to declare the usage of rules from external modules (see :ref:`snakefiles-modules`).
Consider the following example:

.. code-block:: python

    rule a:
        output:
            "test.out"
        shell:
            "echo test > {output}"


    use rule a as b with:
        output:
            "test2.out"


As can be seen, we first declare a rule a, and then we reuse the rule a as rule b, while changing only the output file and keeping everything else the same.
In reality, one will often change more.
Analogously to the ``use rule`` from external modules, any properties of the rule (``input``, ``output``, ``log``, ``params``, ``benchmark``, ``threads``, ``resources``, ``pathvars``, etc.) can be modified, except the actual execution step (``shell``, ``notebook``, ``script``, ``cwl``, or ``run``).
All unmodified properties are inherited from the parent rule.

:ref:`Pathvars <snakefiles-pathvars>` become particularly powerful in combination with such rule inheritance, as they allow to introduce generic items in the parent rule that can be specified out in the child rule:

.. code-block:: python

    rule transform_something:
        input:
            "<results>/<instep>/<per>.txt"
        output:
            "<results>/<outstep>/<per>.txt"
        shell:
            "somecommand {input} > {output}"

    use rule transform_something as something1_to_something2 with:
        pathvars:
            instep="something1",
            outstep="something2",
            per="{sample}"

    use rule transform_something_else as something5_to_something6 with:
        pathvars:
            instep="something5",
            outstep="something6",
            per="{sample}.{replicate}"

In other words, here we define a potentially complex rule only once, and explicitly use it in two different parts of the workflow, even with different kinds of wildcards, all by just configuring the pathvars.


.. important::
    A rule cannot be redefined without renaming it using the ``as`` clause.
    Otherwise, you will have two versions of the same rule, which might be unintended (a common symptom of such unintended repeated uses would be ambiguous rule exceptions thrown by Snakemake).
    However, it is allowed to create **multiple modified versions** of the same rule, as long as each has a **unique name**.
    The only exception is when a rule was previously imported via a general ``use rule * from`` statement, such rules may be **further modified once** under the same final name for convenience (see :ref:`snakefiles-modules`).

.. note::
    Modification of `params` allows the replacement of single keyword arguments.
    Keyword `params` arguments of the original rule that are not defined after `with` are inherited.
    Positional `params` arguments of the original rule are overwritten, if positional `params` arguments are given after `with`.
    All other properties (``input``, ``output``, ...) are entirely overwritten with the values specified after `with`.




.. _snakefiles-procedural-rules:

Procedural rule definition
--------------------------

The name is optional and can be left out, creating an anonymous rule. It can also be overridden by setting a rule's ``name`` attribute.

.. code-block:: python

    for tool in ["bcftools", "freebayes"]:
        rule:
            name:
                f"call_variants_{tool}"
            input:
                f"path/to/{tool}/inputfile"
            output:
                f"path/to/{tool}/outputfile"
            shell:
                f"{tool} {{input}} > {{output}}"


.. sidebar:: Note

    Note that any placeholders in the shell command (like ``{input}``) are always evaluated and replaced
    when the corresponding job is executed, even if they are occurring inside a comment.
    To avoid evaluation and replacement, you have to mask the braces by doubling them,
    i.e. ``{{input}}``.
