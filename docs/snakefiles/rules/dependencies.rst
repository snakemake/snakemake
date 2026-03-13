=================
Rule Dependencies
=================


.. _snakefiles-targets:

Target rules
-------------

By default, Snakemake always wants to execute the first rule in the snakefile.
This gives rise to pseudo-rules at the beginning of the file that can be used to define build-targets similar to GNU Make:

.. code-block:: python

    rule all:
        input:
            expand("{dataset}/file.A.txt", dataset=DATASETS)


Here, for each dataset in a python list ``DATASETS`` defined before, the file ``{dataset}/file.A.txt`` is requested.
In this example, Snakemake recognizes automatically that these can be created by multiple applications of the rule ``complex_conversion`` shown above.

It is possible to overwrite this behavior to use the first rule as a default target, by explicitly marking a rule as being the default target via the ``default_target`` directive:

.. code-block:: python

    rule xy:
        input:
            expand("{dataset}/file.A.txt", dataset=DATASETS)
        default_target: True

Regardless of where this rule appears in the Snakefile, it will be the default target.
Usually, it is still recommended to keep the default target rule (and in fact all other rules that could act as optional targets) at the top of the file, such that it can be easily found.
The ``default_target`` directive becomes particularly useful when :ref:`combining several pre-existing workflows <use_with_modules>`.


Rule dependencies
-----------------

From version 2.4.8 on, rules can also refer to the output of other rules in the Snakefile, e.g.:

.. code-block:: python

    rule a:
        input:  "path/to/input"
        output: "path/to/output"
        shell:  ...

    rule b:
        input:  rules.a.output
        output: "path/to/output/of/b"
        shell:  ...

Importantly, be aware that referring to rule ``a`` here requires that rule ``a`` was defined above rule ``b`` in the file, since the object has to be known already.
This feature also allows us to resolve dependencies that are ambiguous when using filenames.

Note that when the rule you refer to defines multiple output files but you want to require only a subset of those as input for another rule, you should name the output files and refer to them specifically:

.. code-block:: python

    rule a:
        input:  "path/to/input"
        output: a = "path/to/output", b = "path/to/output2"
        shell:  ...

    rule b:
        input:  rules.a.output.a
        output: "path/to/output/of/b"
        shell:  ...


.. _snakefiles-ambiguous-rules:

Handling Ambiguous Rules
------------------------

When two rules can produce the same output file, snakemake cannot decide which one to use without additional guidance. Hence an ``AmbiguousRuleException`` is thrown.
Note: ruleorder is not intended to bring rules in the correct execution order (this is solely guided by the names of input and output files you use), it only helps snakemake to decide which rule to use when multiple ones can create the same output file!
To deal with such ambiguity, provide a ``ruleorder`` for the conflicting rules, e.g.

.. code-block:: python

    ruleorder: rule1 > rule2 > rule3

Here, ``rule1`` is preferred over ``rule2`` and ``rule3``, and ``rule2`` is preferred over ``rule3``.
Only if rule1 and rule2 cannot be applied (e.g. due to missing input files), rule3 is used to produce the desired output file.

Alternatively, rule dependencies (see above) can also resolve ambiguities.

Another (quick and dirty) possibility is to tell snakemake to allow ambiguity via a command line option

.. code-block:: console

    $ snakemake --allow-ambiguity

such that similar to GNU Make always the first matching rule is used. Here, a warning that summarizes the decision of snakemake is provided at the terminal.

.. _snakefiles-scattergather:

Defining scatter-gather processes
---------------------------------

Via Snakemake's powerful and arbitrary Python based aggregation abilities (via the ``expand`` function and arbitrary Python code, see :ref:`here <snakefiles_aggregation>`), scatter-gather workflows are well supported.
Nevertheless, it can sometimes be handy to use Snakemake's specific scatter-gather support, which allows to avoid boilerplate and offers additional configuration options.
Scatter-gather processes can be defined via a global ``scattergather`` directive:

.. code-block:: python

    scattergather:
        split=8

Each process thereby defines a name (here e.g. ``split``) and a default number of scatter items.
Then, scattering and gathering can be implemented by using globally available ``scatter`` and ``gather`` objects:

.. code-block:: python


    rule all:
        input:
            "gathered/all.txt"


    rule split:
        output:
            scatter.split("split/{scatteritem}.txt")
        shell:
            "touch {output}"


    rule intermediate:
        input:
            "split/{scatteritem}.txt"
        output:
            "split/{scatteritem}.post.txt"
        shell:
            "cp {input} {output}"


    rule gather:
        input:
            gather.split("split/{scatteritem}.post.txt")
        output:
            "gathered/all.txt"
        shell:
            "cat {input} > {output}"

Thereby, ``scatter.split("split/{scatteritem}.txt")`` yields a list of paths ``"split/1-of-n.txt"``, ``"split/2-of-n.txt"``, ..., depending on the number ``n`` of scatter items defined.
Analogously, ``gather.split("split/{scatteritem}.post.txt")``, yields a list of paths ``"split/0.post.txt"``, ``"split/1.post.txt"``, ..., which request the application of the rule ``intermediate`` to each scatter item.

The default number of scatter items can be overwritten via the command line interface.
For example

.. code-block:: bash

    snakemake --set-scatter split=2

would set the number of scatter items for the split process defined above to 2 instead of 8.
This allows to adapt parallelization according to the needs of the underlying computing platform and the analysis at hand.

For more complex workflows it's possible to define multiple processes, for example:

.. code-block:: python

    scattergather:
        split_a=8,
        split_b=3,

The calls to ``scatter`` and ``gather`` would need to reference the appropriate process name, e.g. ``scatter.split_a`` and ``gather.split_a`` to use the ``split_a`` settings.


.. _snakefiles-checkpoints:

Data-dependent conditional execution
------------------------------------

From Snakemake 5.4 on, conditional reevaluation of the DAG of jobs based on the content outputs is possible.
The key idea is that rules can be declared as checkpoints, e.g.,

.. code-block:: python

    checkpoint somestep:
        input:
            "samples/{sample}.txt"
        output:
            "somestep/{sample}.txt"
        shell:
            "somecommand {input} > {output}"

Snakemake allows to re-evaluate the DAG after the successful execution of every job spawned from a checkpoint.
For this, every checkpoint is registered by its name in a globally available ``checkpoints`` object.
The ``checkpoints`` object can be accessed by :ref:`input functions <snakefiles-input_functions>`.
Assuming that the checkpoint is named ``somestep`` as above, the output files for a particular job can be retrieved with

.. code-block:: python

    checkpoints.somestep.get(sample="a").output

.. note::

    Note that output files of checkpoints that are accessed via this mechanism will not be marked as temporary.
    Even you try to mark them as temporary, Snakemake will ignore the label and keep the output files of the checkpoint.
    Reruns will not be triggered if the output file do not exist.

Thereby, the ``get`` method throws ``snakemake.exceptions.IncompleteCheckpointException`` if the checkpoint has not yet been executed for these particular wildcard value(s).
Inside an input function, the exception will be automatically handled by Snakemake, and leads to a re-evaluation after the checkpoint has been successfully passed.

To illustrate the possibilities of this mechanism, consider the following complete example:

.. code-block:: python

    # a target rule to define the desired final output
    rule all:
        input:
            "aggregated/sample1.txt",
            "aggregated/sample2.txt"


    # generate per-sample input files; filenames are sample IDs,
    # while file content ("a" or "b") controls downstream branching
    rule generate_sample_input:
        output:
            "samples/{sample}.txt"
        run:
            import random

            with open(output[0], "w") as f:
                f.write(random.choice(["a", "b"]))


    # the checkpoint that shall trigger re-evaluation of the DAG
    checkpoint somestep:
        input:
            "samples/{sample}.txt"
        output:
            "somestep/{sample}.txt"
        shell:
            # propagate generated value into checkpoint output
            "cp {input} {output}"


    # intermediate rule
    rule intermediate:
        input:
            "somestep/{sample}.txt"
        output:
            "post/{sample}.txt"
        shell:
            "touch {output}"


    # alternative intermediate rule
    rule alt_intermediate:
        input:
            "somestep/{sample}.txt"
        output:
            "alt/{sample}.txt"
        shell:
            "touch {output}"


    # input function for the rule aggregate
    def aggregate_input(wildcards):
        # decision based on content of output file
        # Important: use the method open() of the returned file!
        # This way, Snakemake is able to automatically download the file if it is generated in
        # a cloud environment without a shared filesystem.
        with checkpoints.somestep.get(sample=wildcards.sample).output[0].open() as f:
            if f.read().strip() == "a":
                return "post/{sample}.txt"
            else:
                return "alt/{sample}.txt"


    rule aggregate:
        input:
            aggregate_input
        output:
            "aggregated/{sample}.txt"
        shell:
            "touch {output}"

As can be seen, the rule aggregate uses an input function.

.. note::

    You don't need to use the checkpoint mechanism to determine parameter or resource values of downstream rules that would be based on the output of previous rules.
    In fact, it won't even work because the checkpoint mechanism is only considered for input functions.
    Instead, you can simply use normal parameter or resource functions that just assume that those output files are there. Snakemake will evaluate them immediately before
    the job is scheduled, when the required files from upstream rules are already present.

Inside the function, we first retrieve the output files of the checkpoint ``somestep`` with the wildcards, passing through the value of the wildcard sample.
Upon execution, if the checkpoint is not yet complete, Snakemake will record ``somestep`` as a direct dependency of the rule ``aggregate``.
Once ``somestep`` has finished for a given sample, the input function will automatically be re-evaluated and the method ``get`` will no longer return an exception.
Instead, the output file will be opened, and depending on its contents either ``"post/{sample}.txt"`` or ``"alt/{sample}.txt"`` will be returned by the input function.
This way, the DAG becomes conditional on some produced data.

It is also possible to use checkpoints for cases where the output files are unknown before execution.
Consider the following example where an arbitrary number of files is generated by a rule before being aggregated:

.. code-block:: python

    # a target rule to define the desired final output
    rule all:
        input:
            "aggregated.txt"


    # the checkpoint that shall trigger re-evaluation of the DAG
    # an number of file is created in a defined directory
    checkpoint somestep:
        output:
            directory("my_directory/")
        shell:'''
        mkdir my_directory/
        cd my_directory
        for i in 1 2 3; do touch $i.txt; done
        '''



    # input function for rule aggregate, return paths to all files produced by the checkpoint 'somestep'
    def aggregate_input(wildcards):
        checkpoint_output = checkpoints.somestep.get(**wildcards).output[0]
        return expand("my_directory/{i}.txt",
                    i=glob_wildcards(os.path.join(checkpoint_output, "{i}.txt")).i)


    rule aggregate:
        input:
            aggregate_input
        output:
            "aggregated.txt"
        shell:
            "cat {input} > {output}"

Because the number of output files is unknown beforehand, the checkpoint only defines an output :ref:`directory <snakefiles-directory_output>`.
This time, instead of explicitly writing

.. code-block:: python

    checkpoints.somestep.get(sample=wildcards.sample).output[0]

we use the shorthand

.. code-block:: python

    checkpoints.somestep.get(**wildcards).output[0]

which automatically unpacks the wildcards as keyword arguments (this is standard python argument unpacking).
If the checkpoint has not yet been executed, accessing ``checkpoints.somestep.get(**wildcards)`` ensures that Snakemake records the checkpoint as a direct dependency of the rule ``aggregate``.
Upon completion of the checkpoint, the input function is re-evaluated, and the code beyond its first line is executed.
Here, we retrieve the values of the wildcard ``i`` based on all files named ``{i}.txt`` in the output directory of the checkpoint.
Because the wildcard ``i`` is evaluated only after completion of the checkpoint, it is necessary to use ``directory`` to declare its output, instead of using the full wildcard patterns as output.

A more practical example building on the previous one is a clustering process with an unknown number of clusters for different samples, where each cluster shall be saved into a separate file.
In this example the clusters are being processed by an intermediate rule before being aggregated:

.. code-block:: python

    # a target rule to define the desired final output
    rule all:
        input:
            "aggregated/a.txt",
            "aggregated/b.txt"


    # the checkpoint that shall trigger re-evaluation of the DAG
    checkpoint clustering:
        input:
            "samples/{sample}.txt"
        output:
            clusters=directory("clustering/{sample}")
        shell:
            "mkdir clustering/{wildcards.sample}; "
            "for i in 1 2 3; do echo $i > clustering/{wildcards.sample}/$i.txt; done"


    # an intermediate rule
    rule intermediate:
        input:
            "clustering/{sample}/{i}.txt"
        output:
            "post/{sample}/{i}.txt"
        shell:
            "cp {input} {output}"


    def aggregate_input(wildcards):
        checkpoint_output = checkpoints.clustering.get(**wildcards).output[0]
        return expand("post/{sample}/{i}.txt",
                sample=wildcards.sample,
                i=glob_wildcards(os.path.join(checkpoint_output, "{i}.txt")).i)


    # an aggregation over all produced clusters
    rule aggregate:
        input:
            aggregate_input
        output:
            "aggregated/{sample}.txt"
        shell:
            "cat {input} > {output}"

Here a new directory will be created for each sample by the checkpoint.
After completion of the checkpoint, the ``aggregate_input`` function is re-evaluated as previously.
The values of the wildcard ``i`` is this time used to expand the pattern ``"post/{sample}/{i}.txt"``, such that the rule ``intermediate`` is executed for each of the determined clusters.
