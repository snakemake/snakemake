=====================
Input/Output Handling
=====================


.. _snakefiles_aggregation:

Aggregation
-----------

Input files can be Python lists, allowing to easily aggregate over parameters or samples:

.. code-block:: python

    rule aggregate:
        input:
            ["{dataset}/a.txt".format(dataset=dataset) for dataset in DATASETS]
        output:
            "aggregated.txt"
        shell:
            ...

While the above expression can be very powerful as arbitrary Python code can be used, Snakemake offers
various helper functions to simplify aggregations (see :ref:`snakefiles-input_helpers`).

.. _snakefiles-input_functions:

Input functions
---------------

Instead of specifying strings or lists of strings as input files, snakemake can also make use of functions that return single **or** lists of input files:

.. code-block:: python

    def myfunc(wildcards):
        return [... a list of input files depending on given wildcards ...]

    rule:
        input:
            myfunc
        output:
            "someoutput.{somewildcard}.txt"
        shell:
            "..."

The function has to accept a single argument that will be the wildcards object generated from the application of the rule to create some requested output files.
Note that you can also use `lambda expressions <https://docs.python.org/3/tutorial/controlflow.html#lambda-expressions>`_ instead of full function definitions.
By this, rules can have entirely different input files (both in form and number) depending on the inferred wildcards. E.g. you can assign input files that appear in entirely different parts of your filesystem based on some wildcard value and a dictionary that maps the wildcard value to file paths.

.. note::

    Input functions can themselves return input functions again (this also holds for functions given to params and resources.)
    Such nested evaluation is allowed for a depth up to 10. Afterwards, an exception will be thrown.

In addition to a single wildcards argument, input functions can optionally take a ``groupid`` (with exactly that name) as second argument, see :ref:`snakefiles_group-local` for details.

Finally, when implementing the input function, it is best practice to make sure that it can properly handle all possible wildcard values your rule can have.
In particular, input files should not be combined with very general rules that can be applied to create almost any file: Snakemake will try to apply the rule, and will report the exceptions of your input function as errors.

For a practical example, see the :ref:`tutorial` (:ref:`tutorial-input_functions`).

.. _snakefiles-unpack:

Input Functions and ``unpack()``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In some cases, you might want to have your input functions return named input files.
This can be done by having them return ``dict()`` objects with the names as the dict keys and the file names as the dict values and using the ``unpack()`` keyword.

.. code-block:: python

    def myfunc(wildcards):
        return {'foo': '{wildcards.token}.txt'.format(wildcards=wildcards)}

    rule:
        input:
            unpack(myfunc)
        output:
            "someoutput.{token}.txt"
        shell:
            "..."

Note that ``unpack()`` is only necessary for input functions returning ``dict``.
While it also works for ``list``, remember that lists (and nested lists) of strings are automatically flattened.

Also note that if you do not pass in a *function* into the input list but you directly *call a function* then you shouldn't use ``unpack()``.
Here, you can simply use Python's double-star (``**``) operator for unpacking the parameters.

Note that as Snakefiles are translated into Python for execution, the same rules as for using the `star and double-star unpacking Python operators <https://docs.python.org/3/tutorial/controlflow.html#unpacking-argument-lists>`_ apply.
These restrictions do not apply when using ``unpack()``.

.. code-block:: python

    def myfunc1():
        return ['foo.txt']

    def myfunc2():
        return {'foo': 'nowildcards.txt'}

    rule:
        input:
            *myfunc1(),
            **myfunc2(),
        output:
            "..."
        shell:
            "..."


.. _snakefiles_protected_temp:

Protected and Temporary Files
-----------------------------

A particular output file may require a huge amount of computation time. Hence one might want to protect it against accidental deletion or overwriting. Snakemake allows this by marking such a file as ``protected``:

.. code-block:: python

    rule NAME:
        input:
            "path/to/inputfile"
        output:
            protected("path/to/outputfile")
        shell:
            "somecommand {input} {output}"

A protected file will be write-protected after the rule that produces it is completed.

Further, an output file marked as ``temp`` is deleted after all rules that use it as an input are completed:

.. code-block:: python

    rule NAME:
        input:
            "path/to/inputfile"
        output:
            temp("path/to/outputfile")
        shell:
            "somecommand {input} {output}"

.. _snakefiles-directory_output:

Auto-grouping ``temp`` Files in Remote Execution
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For performance reasons, it is sometimes useful to write intermediate files on a faster storage, e.g., attached locally on the cluster compute node rather than shared over the network (and thus neither visible to the main snakemake process that submits jobs to the cluster, nor to other nodes of the cluster).
Snakemake (since version 9.0) allows files marked as ``temp`` to use the option ``group_jobs`` to indicate that rules creating and consuming them should be automatically :ref:`grouped  <job_grouping>` together so Snakemake will schedule them to run on the same physical node:

.. code-block:: python

    rule NAME1:
        input:
            "path/to/inputfile"
        output:
            temp("path/to/intermediatefile", group_jobs=True)
        shell:
            "somecommand {input} {output}"

    rule NAME2:
        input:
            "path/to/intermediatefile"
        output:
            "path/to/outputfile"
        shell:
            "someothercommand {input} {output}"

Directories as Outputs
----------------------

Sometimes it can be convenient to have directories, rather than files, as outputs of a rule.
As of version 5.2.0, directories as outputs have to be explicitly marked with ``directory``.
This is primarily for safety reasons; since all outputs are deleted before a job is executed, we don't want to risk deleting important directories if the user makes some mistake.
Marking the output as ``directory`` makes the intent clear, and the output can be safely removed.
Another reason comes down to how modification time for directories work.
The modification time on a directory changes when a file or a subdirectory is added, removed or renamed.
This can easily happen in not-quite-intended ways, such as when Apple macOS or MS Windows add ``.DS_Store`` or ``thumbs.db`` files to store parameters for how the directory contents should be displayed.
When the ``directory`` flag is used a hidden file called ``.snakemake_timestamp`` is created in the output directory, and the modification time of that file is used when determining whether the rule output is up to date or if it needs to be rerun.
Always consider if you can't formulate your workflow using normal files before resorting to using ``directory()``.

.. code-block:: python

    rule NAME:
        input:
            "path/to/inputfile"
        output:
            directory("path/to/outputdir")
        shell:
            "somecommand {input} {output}"

.. sidebar:: Caution on Directory Outputs

    Note that because a directory marked as the output of a job will be deleted before that job executes (in order to start with a clean state), other jobs should not create output files within that directory, as these may be inadvertently deleted as well.

Ignoring Timestamps
-------------------

For determining whether output files have to be re-created, Snakemake checks whether the file modification date (i.e. the timestamp) of any input file of the same job is newer than the timestamp of the output file.
This behavior can be overridden by marking an input file as ``ancient``.
The timestamp of such files is ignored and always assumed to be older than any of the output files:

.. code-block:: python

    rule NAME:
        input:
            ancient("path/to/inputfile")
        output:
            "path/to/outputfile"
        shell:
            "somecommand {input} {output}"

Here, this means that the file ``path/to/outputfile`` will not be triggered for re-creation after it has been generated once, even when the input file is modified in the future.
Note that any flag that forces re-creation of files still also applies to files marked as ``ancient``.

.. _snakefiles_ensure:

Ensuring Output File Properties
-------------------------------

It is possible to annotate certain additional criteria for output files to be ensured after they have been generated successfully.
For example, this can be used to check for output files to be non-empty, or to compare them against a given sha256 checksum.
If this functionality is used, Snakemake will check such annotated files before considering a job to be successful.
Non-emptyness can be checked as follows:

.. code-block:: python

    rule NAME:
        output:
            ensure("test.txt", non_empty=True)
        shell:
            "somecommand {output}"

Above, the output file ``test.txt`` is marked as non-empty.
If the command ``somecommand`` happens to generate an empty output,
the job will fail with an error listing the unexpected empty file.

A sha256 (or md5 or sha1) checksum can be compared as follows (using corresponding keyword arguments ``sha256=``, ``md5=``, or ``sha1=``).:

.. code-block:: python

    my_checksum = "9f86d081884c7d659a2feaa0c55ad015a3bf4f1b2b0b822cd15d6c15b0f00a08"

    rule NAME:
        output:
            ensure("test.txt", sha256=my_checksum)
        shell:
            "somecommand {output}"

In addition to providing the checksum as plain string, it is possible to provide a pointer to a function (similar to :ref:`input functions <snakefiles-input_functions>`).
The function has to accept a single argument that will be the wildcards object generated from the application of the rule to create some requested output files:

.. code-block:: python

    def get_checksum(wildcards):
        # e.g., look up the checksum with the value of the wildcard sample
        # in some dictionary
        return my_checksums[wildcards.sample]

    rule NAME:
        output:
            ensure("test/{sample}.txt", sha256=get_checksum)
        shell:
            "somecommand {output}"


Note that you can also use `lambda expressions <https://docs.python.org/3/tutorial/controlflow.html#lambda-expressions>`_ instead of full function definitions.

Often, it is a good idea to combine ``ensure`` annotations with :ref:`retry definitions <snakefiles_retries>`, e.g. for retrying upon invalid checksums or empty files.

Shadow Rules
------------

Shadow rules result in each execution of the rule to be run in isolated temporary directories.
This "shadow" directory contains symlinks to files and directories in the current workdir.
This is useful for running programs that generate lots of unused files which you don't want to manually cleanup in your snakemake workflow.
It can also be useful if you want to keep your workdir clean while the program executes,
or simplify your workflow by not having to worry about unique filenames for all outputs of all rules.

By setting ``shadow: "shallow"``, the top level files and directories are symlinked,
so that any relative paths in a subdirectory will be real paths in the filesystem.
The setting ``shadow: "full"`` fully shadows the entire subdirectory structure of the current workdir.
The setting ``shadow: "minimal"`` only symlinks the inputs to the rule,
and ``shadow: "copy-minimal"`` copies the inputs instead of just creating symlinks.
Once the rule successfully executes, the output file will be moved if necessary to the real path as indicated by ``output``.

Typically, you will not need to modify your rule for compatibility with ``shadow``,
unless you reference parent directories relative to your workdir in a rule.

.. code-block:: python

    rule NAME:
        input: "path/to/inputfile"
        output: "path/to/outputfile"
        shadow: "shallow"
        shell: "somecommand --other_outputs other.txt {input} {output}"

Shadow directories are stored one per rule execution in ``.snakemake/shadow/``,
and are cleared on successful execution.
Consider running with the ``--cleanup-shadow`` argument every now and then
to remove any remaining shadow directories from aborted jobs.
The base shadow directory can be changed with the ``--shadow-prefix`` command line argument.



.. _snakefiles-piped-output:

Piped output
------------

From Snakemake 5.0 on, it is possible to mark output files as pipes, via the ``pipe`` flag, e.g.:

.. code-block:: python

  rule all:
      input:
          expand("test.{i}.out", i=range(2))


  rule a:
      output:
          pipe("test.{i}.txt")
      shell:
          "for i in {{0..2}}; do echo {wildcards.i} >> {output}; done"


  rule b:
      input:
          "test.{i}.txt"
      output:
          "test.{i}.out"
      shell:
          "grep {wildcards.i} < {input} > {output}"

If an output file is marked to be a pipe, then Snakemake will first create a `named pipe <https://en.wikipedia.org/wiki/Named_pipe>`_ with the given name and then execute the creating job simultaneously with the consuming job, inside a **group job** (see above).
This works in all execution modes, local, cluster, and cloud.
Naturally, a pipe output may only have a single consumer.
It is possible to combine explicit group definition as above with pipe outputs.
Thereby, pipe jobs can live within, or (automatically) extend existing groups.
However, the two jobs connected by a pipe may not exist in conflicting groups.

As with other groups, Snakemake will automatically calculate the required resources for the group job (see :ref:`resources <snakefiles-resources>`.



.. _snakefiles_continuous_input:

Continuously updated input
--------------------------

From Snakemake 8.2 on, it is possible to define rules that continuously accept new input files during workflow execution.
This is useful for scenarios like streaming data analysis.
The feature works by defining a synchronized Python queue for obtaining input files via the helper function ``from_queue``:

.. code-block:: python

    rule myrule:
        input:
            from_queue(all_results, finish_sentinel=...)
        ...

Rules with input marked as ``from_queue`` may not define any wildcards.
When new items arrive in the queue:

1. The input files list for the rule is updated
2. The DAG of jobs is updated, potentially generating new dependencies for the rule
3. Any dependent rules that need to process the new input files are automatically created and executed

It is required to define a finish sentinel, which is a special value that signals the end of the queue.
Once the finish sentinel is encountered, Snakemake will allow all remaining dependent jobs to finish and complete execution of the workflow.

Consider the following complete toy example:

.. code-block:: python

    import threading, queue, time

    # the finish sentinel
    finish_sentinel = object()
    # a synchronized queue for the input files
    all_results = queue.Queue()

    # a thread that fills the queue with input files to be considered
    def update_results():
        try:
            for i in range(10):
                all_results.put(f"test{i}.txt")
                time.sleep(1)
            all_results.put(finish_sentinel)
            all_results.join()
        except (KeyboardInterrupt, SystemExit):
            return

    update_thread = threading.Thread(target=update_results)
    update_thread.start()


    # target rule which will be continuously updated until the queue is finished
    rule all:
        input:
            from_queue(all_results, finish_sentinel=finish_sentinel)


    # job that generates the requested input files
    rule generate:
        output:
            "test{i}.txt"
        shell:
            "echo {wildcards.i} > {output}"

.. _snakefiles_update_output:

Updating Existing Output Files
------------------------------

By default, Snakemake deletes already existing output files before a job is executed.
This is usually very convenient, because many tools will fail if their output files already exist.
However, from Snakemake 8.7 on, it is possible to declare an output file/directory to be updated by a job instead of rewritten from scratch.
Consider the following example:

.. code-block:: python

    rule update:
        input:
            "in.txt"
        output:
            update("test.txt")
        shell:
            "echo test >> {output}"


Here, the statement ``test`` is appended to the output file ``test.txt``.
Hence, we declare it as being updated via the ``update`` flag.
This way, Snakemake will not delete the file or directory before the job is executed.
Furthermore, Snakemake will restore the previous version of the file/directory if the job fails.

If such a file/directory has to be considered as input **before the update** for another rule
it can be marked as ``before_update``.
This ensures that Snakemake does not search for a producing job but instead considers the file as is on disk or in the storage:

.. code-block:: python

    rule do_something:
        input:
            before_update("test.txt")
        output:
            "in.txt"
        shell:
            "cp {input} {output}"

    rule update:
        input:
            "in.txt"
        output:
            update("test.txt")
        shell:
            "echo test >> {output}"

As can be seen, this way it is even possible to break a cyclic dependency.
An important helper for setting up the logic of ``before_update`` is the :ref:`exists function <snakefiles-semantic-helpers-exists>`, which allows to e.g. condition the consideration of the file that shall be used before the update by its actual existence before the update.



Flag Files
----------

Sometimes it is necessary to enforce some rule execution order without real file dependencies. This can be achieved by "touching" empty files that denote that a certain task was completed. Snakemake supports this via the `touch` flag:

.. code-block:: python

    rule all:
        input: "mytask.done"

    rule mytask:
        output: touch("mytask.done")
        shell: "mycommand ..."

With the ``touch`` flag, Snakemake touches (i.e. creates or updates) the file ``mytask.done`` after ``mycommand`` has finished successfully.


.. _snakefiles-aux_source_files:

Accessing auxiliary source files
--------------------------------

Snakemake workflows can refer to various other source files via paths relative to the current Snakefile.
This happens for example with the :ref:`script directive <snakefiles-external_scripts>` or the :ref:`conda directive <integrated_package_management>`.
Sometimes, it is necessary to access further source files that are in a directory relative to the current Snakefile.
Since workflows can be imported from remote locations (e.g. when using :ref:`modules <snakefiles-modules>`), it is important to not do this manually, so that Snakemake has the chance to cache these files locally before they are accessed.
This can be achieved by accessing their path via the ``workflow.source_path``, which (a) computes the correct path relative to the current Snakefile such that the file can be accessed from any working directory, and (b) downloads remote files to a local cache:

.. code-block:: python

    rule a:
        input:
            json=workflow.source_path("../resources/test.json")
        output:
            "test.out"
        shell:
            "somecommand {input.json} > {output}"

.. note::
    Note that if such source paths are specified as input files, they are automatically considered to be non-storage files.
    This means that Snakemake will not try to map them to an eventually specified default storage provider (see :ref:`storage-support`).
    Further, note that ``workflow.source_path`` should not be used from ``params:`` but only from ``input:``. The reason is that it returns a cached path that may change between Snakemake runs, thereby triggering spurious reruns if referred via ``params:`` (since Snakemake would think that the parameter has changed.


.. _snakefiles_default_flags:

Setting default flags
---------------------

Snakemake allows the annotation of input and output files via so-called flags (see e.g. :ref:`snakefiles_protected_temp`).
Sometimes, it can be useful to define that a certain flag shall be applied to all input or output files of a workflow.
This can be achieved via the global ``inputflags`` and ``outputflags`` directives.
Consider the following example:

.. code-block:: python

    outputflags:
        temp

    rule a:
        output:
            "test.out"
        shell:
            "echo test > {output}"

Would automatically mark the output file of rule ``a`` as temporary.
The most convenient use case of this mechanism occurs in combination with :ref:`access pattern annotation <storage-access-patterns>`.
In this case, the default access pattern can be set globally for all output files of a workflow.
Only a few cases that differ have then to deal with explicit access pattern annotation (see :ref:`storage-access-patterns` for an example).
Whenever a rule defines a flag for a file, this flag will override the default flag of the same kind or any contradicting default flags (e.g. ``temp`` will override ``protected``).

Such default input and output flag specifications are always valid for all rules that follow them in the workflow definition.
Importantly, they are also "namespaced" per module, meaning that ``inputflags`` and ``outputflags`` directives in a module only apply to the rules defined in that module.
