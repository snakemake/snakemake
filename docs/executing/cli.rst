.. _executable:

======================
Command line interface
======================

This part of the documentation describes the ``snakemake`` executable.  Snakemake
is primarily a command-line tool, so the ``snakemake`` executable is the primary way
to execute, debug, and visualize workflows.

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
Further, the reason for each rule execution can be printed via


.. code-block:: console

    $ snakemake -n -r

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
This can be particularly handy when used in combination with :ref:`cluster execution <cluster>`.

Dealing with very large workflows
---------------------------------

If your workflow has a lot of jobs, Snakemake might need some time to infer the dependencies (the job DAG) and which jobs are actually required to run.
The major bottleneck involved is the filesystem, which has to be queried for existence and modification dates of files.
To overcome this issue, Snakemake allows to run large workflows in batches.
This way, fewer files have to be evaluated at once, and therefore the job DAG can be inferred faster.
By running

.. code-block:: console

    $ snakemake --cores 4 --batch myrule=1/3

you instruct to only compute the first of three batches of the inputs of the rule `myrule`.
To generate the second batch, run

.. code-block:: console

    $ snakemake --cores 4 --batch myrule=2/3

Finally, when running


.. code-block:: console

    $ snakemake --cores 4 --batch myrule=3/3

Snakemake will process beyond the rule `myrule`, because all of its input files have been generated, and complete the workflow.
Obviously, a good choice of the rule to perform the batching is a rule that has a lot of input files and upstream jobs, for example a central aggregation step within your workflow.
We advice all workflow developers to inform potential users of the best suited batching rule.

.. _profiles:

--------
Profiles
--------

Adapting Snakemake to a particular environment can entail many flags and options.
Therefore, since Snakemake 4.1, it is possible to specify a configuration profile
to be used to obtain default options:

.. code-block:: console

   $ snakemake --profile myprofile

Here, a folder ``myprofile`` is searched in per-user and global configuration directories (on Linux, this will be ``$HOME/.config/snakemake`` and ``/etc/xdg/snakemake``, you can find the answer for your system via ``snakemake --help``).
Alternatively, an absolute or relative path to the folder can be given.
The profile folder is expected to contain a file ``config.yaml`` that defines default values for the Snakemake command line arguments.
For example, the file

.. code-block:: yaml

    cluster: qsub
    jobs: 100

would setup Snakemake to always submit to the cluster via the ``qsub`` command, and never use more than 100 parallel jobs in total.
The profile can be used to set a default for each option of the Snakemake command line interface.
For this, option ``--someoption`` becomes ``someoption: `` in the profile.
If options accept multiple arguments these must be given as YAML list in the profile.
Under https://github.com/snakemake-profiles/doc, you can find publicly available profiles.
Feel free to contribute your own.

The profile folder can additionally contain auxilliary files, e.g., jobscripts, or any kind of wrappers.
See https://github.com/snakemake-profiles/doc for examples.


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
   :module: snakemake
   :func: get_argument_parser
   :prog: snakemake

   All command line options can be printed by calling ``snakemake -h``.

.. _getting_started-bash_completion:

---------------
Bash Completion
---------------

Snakemake supports bash completion for filenames, rulenames and arguments.
To enable it globally, just append

.. code-block:: bash

    `snakemake --bash-completion`

including the accents to your ``.bashrc``.
This only works if the ``snakemake`` command is in your path.
