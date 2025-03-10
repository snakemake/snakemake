.. _executable:

======================
Command line interface
======================

This part of the documentation describes the ``snakemake`` executable.  Snakemake
is primarily a command-line tool, so the ``snakemake`` executable is the primary way
to execute, debug, and visualize workflows.

.. user_manual-snakemake_envvars:

-------------------------------
Important environment variables
-------------------------------

Snakemake caches source files for performance and reproducibility.
The location of this cache is determined by the `appdirs <https://github.com/ActiveState/appdirs>`_ package.
If you want to change the location on a unix/linux system, you can define an override path via the environment variable ``XDG_CACHE_HOME``.

.. user_manual-snakemake_options:

-----------------------------
Useful Command Line Arguments
-----------------------------

If called with the number of cores to use, i.e.

.. code-block:: console

    $ snakemake --cores 1

Snakemake tries to execute the workflow specified in a file called ``Snakefile`` in the same directory (the Snakefile can be given via the parameter ``-s``).

By issuing

.. code-block:: console

    $ snakemake -n

a dry-run can be performed.
This is useful to test if the workflow is defined properly and to estimate the amount of needed computation.

Importantly, Snakemake can automatically determine which parts of the workflow can be run in parallel.
By specifying more than one available core, i.e.

.. code-block:: console

    $ snakemake --cores 4

one can tell Snakemake to use up to 4 cores and solve a binary knapsack problem to optimize the scheduling of jobs.
If the number is omitted (i.e., only ``--cores`` is given), the number of used cores is determined as the number of available CPU cores in the machine.

Snakemake workflows usually define the number of used threads of certain rules. Sometimes, it makes sense to overwrite the defaults given in the workflow definition.
This can be done by using the ``--set-threads`` argument, e.g.,

.. code-block:: console

    $ snakemake --cores 4 --set-threads myrule=2

would overwrite whatever number of threads has been defined for the rule ``myrule`` and use ``2`` instead.
Similarly, it is possible to overwrite other resource definitions in rules, via

.. code-block:: console

    $ snakemake --cores 4 --set-resources myrule:partition="foo"

Both mechanisms can be particularly handy when used in combination with :ref:`non-local execution <non-local-exec>`.

.. _non-local-exec:

Non-local execution
^^^^^^^^^^^^^^^^^^^

Non-local execution on cluster or cloud infrastructure is implemented via plugins.
The `Snakemake plugin catalog <https://snakemake.github.io/snakemake-plugin-catalog>`_ lists available plugins and their documentation.

Dealing with very large workflows
---------------------------------

If your workflow has a lot of jobs, Snakemake might need some time to infer the dependencies (the job DAG) and which jobs are actually required to run.
The major bottleneck involved is the filesystem, which has to be queried for existence and modification dates of files.
To overcome this issue, Snakemake allows to run large workflows in batches.
This way, fewer files have to be evaluated at once, and therefore the job DAG can be inferred faster.
By running

.. code-block:: console

    $ snakemake --cores 4 --batch myrule=1/3

you instruct to only compute the first of three batches of the inputs of the rule ``myrule``.
To generate the second batch, run

.. code-block:: console

    $ snakemake --cores 4 --batch myrule=2/3

Finally, when running


.. code-block:: console

    $ snakemake --cores 4 --batch myrule=3/3

Snakemake will process beyond the rule ``myrule``, because all of its input files have been generated, and complete the workflow.
Obviously, a good choice of the rule to perform the batching is a rule that has a lot of input files and upstream jobs, for example a central aggregation step within your workflow.
We advice all workflow developers to inform potential users of the best suited batching rule.

.. _profiles:

--------
Profiles
--------

Adapting Snakemake to a particular environment can entail many flags and options.
Therefore, since Snakemake 4.1, it is possible to specify configuration profiles
to be used to obtain default options.
Since Snakemake 7.29, two kinds of profiles are supported:

* A **global profile** that is defined in a system-wide or user-specific configuration directory (on Linux, this will be ``$HOME/.config/snakemake`` and ``/etc/xdg/snakemake``, you can find the answer for your system via ``snakemake --help``).
* A **workflow specific profile** (introduced in Snakemake 7.29) that is defined via a flag (``--workflow-profile``) or searched in a default location (``profiles/default``) in the working directory or next to the Snakefile.

The workflow specific profile is meant to be used to define default options for a particular workflow, like providing constraints for certain custom resources the workflow uses (e.g. ``api_calls``) or overwriting the threads and resource definitions of individual rules without modifying the workflow code itself.
In contrast, the global profile is meant to be used to define default options for a particular environment, like the default cluster submission command or the default number of jobs to run in parallel.

For example, the command

.. code-block:: console

   $ snakemake --profile myprofile

would expect a folder ``myprofile`` in per-user and global configuration directories (on Linux, this will be ``$HOME/.config/snakemake`` and ``/etc/xdg/snakemake``, you can find the answer for your system via ``snakemake --help``).
Alternatively, an absolute or relative path to the profile folder can be given.
The default profile to use when no ``--profile`` argument is specified can also be set via the environment variable ``SNAKEMAKE_PROFILE``,
e.g. by specifying ``export SNAKEMAKE_PROFILE=myprofile`` in your ``~/.bashrc`` or the system wide shell defaults means that the ``--profile`` flag can be omitted.
In order unset the profile defined by this environment variable for individual runs without specifying and alternative profile you can provide the special value ``none``, i.e. ``--profile none``.

The profile folder is expected to contain a configuration file that file that defines default values for the Snakemake command line arguments.
The file has to be named ``config.vX+.yaml`` with ``X`` denoting the minimum supported Snakemake major version (e.g. ``config.v8+.yaml``).
As fallback, it is also possible to provide a version agnostic ``config.yaml`` that matches any Snakemake version.
For example, the file

.. code-block:: yaml

    executor: slurm
    jobs: 100

would setup Snakemake to always submit to the SLURM cluster middleware and never use more than 100 parallel jobs in total.
The profile can be used to set a default for each option of the Snakemake command line interface.
For this, option ``--someoption`` becomes ``someoption:`` in the profile.
The profile folder can additionally contain auxiliary files, e.g., jobscripts, or any kind of wrappers. See https://github.com/snakemake-profiles/doc for examples.
If options accept multiple arguments these must be given as YAML list in the profile.
If options expect structured arguments (like ``--default-resources RESOURCE=VALUE``, ``--set-threads RULE=VALUE``, or ``--set-resources RULE:RESOURCE=VALUE``), those can be given as strings in the expected forms, i.e.

.. code-block:: yaml

    default-resources: mem_mb=200
    set-threads: myrule=5
    set-resources: myrule:mem=500MB

or as YAML maps, which is easier to read:

.. code-block:: yaml

    default-resources:
        mem_mb: 200
    set-threads:
        myrule: 5
    set-resources:
        myrule:
            mem: 500MB

All of these resource specifications can also be made dynamic, by using expressions and certain variables that are available.
For details of the variables you can use, refer to the callable signatures given in the
documentation sections on the specification of :ref:`threads <snakefiles-threads>`
and :ref:`dynamic resources <snakefiles-dynamic-resources>`.
These enable ``config.yaml`` entries like:

.. code-block:: yaml

    default-resources:
        mem_mb: max(1.5 * input.size_mb, 100)
    set-threads:
        myrule: max(input.size_mb / 5, 2)
    set-resources:
        myrule:
            mem_mb: attempt * 200


Setting resources or threads via the profile is of course rather a job for the workflow profile instead of the global profile (as such settings are likely workflow specific).

Values in profiles can make use of globally available environment variables, e.g. the ``$USER`` variable.
For example, the following would set the default prefix for storing local copies of remote storage files to a user specific directory

.. code-block:: yaml

    local-storage-prefix: /local/work/$USER/snakemake-scratch

Any such environment variables are automatically expanded when evaluating the profile.

Under https://github.com/snakemake-profiles/doc, you can find publicly available global profiles (e.g. for cluster systems).
Feel free to contribute your own.
Workflow specific profiles are either not shared at all, or can be distributed along with the workflow itself where it makes sense.
For example, when the workflow has its Snakefile at ``workflow/Snakefile``, the profile config should be placed at ``workflow/profiles/default/config.yaml``.


Use templating in profiles
^^^^^^^^^^^^^^^^^^^^^^^^^^

In Snakemake 7.30 or newer, when the profile starts with

.. code-block:: yaml

    __use_yte__: true

It will be treated as a `YTE template <https://yte-template-engine.github.io>`_ and parsed accordingly.
This can be handy to e.g. define values inside of the profile that are based on environment variables.
For example, admins could use this to define user-specific settings.
Another application would be the uniform redefinition of resource requirements for a larger set of rules in a workflow profile (see above).
However, it should be noted that templated profiles are harder to keep free of errors and the profile author has to make sure that they always work correctly for the user.


.. _getting_started-visualization:

-------------
Visualization
-------------

To visualize the workflow, one can use the option ``--dag``.
This creates a representation of the DAG in the graphviz dot language which has to be postprocessed by the graphviz tool ``dot``.
E.g. to visualize the DAG that would be executed, you can issue:

.. code-block:: console

    $ snakemake --dag | dot | display

For saving this to a file, you can specify the desired format:

.. code-block:: console

    $ snakemake --dag | dot -Tpdf > dag.pdf

To visualize the whole DAG regardless of the eventual presence of files, the ``forceall`` option can be used:

.. code-block:: console

    $ snakemake --forceall --dag | dot -Tpdf > dag.pdf

Of course the visual appearance can be modified by providing further command line arguments to ``dot``.

**Note:** The DAG is printed in DOT format straight to the standard output, along with other ``print`` statements you may have in your Snakefile. Make sure to comment these other ``print`` statements so that ``dot`` can build a visual representation of your DAG.


.. _all_options:

-----------
All Options
-----------

.. argparse::
   :module: snakemake.cli
   :func: get_argument_parser
   :prog: snakemake

   All command line options can be printed by calling ``snakemake -h``.

