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
The location of this cache is determined by the `platformdirs <https://platformdirs.readthedocs.io/>`_ package.
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

Non-local execution on cluster or cloud infrastructure is implemented via :ref:`executor plugins <executors>`.
The `Snakemake plugin catalog <https://snakemake.github.io/snakemake-plugin-catalog>`__ lists available plugins and their documentation.
In general, the configuration boils down to specifying an executor plugin (e.g. for SLURM or Kubernetes) and, if needed, a :ref:`storage <default_storage>` plugin (e.g. in order to use S3 for input and output files or in order to efficiently use a shared network filesystem).
For maximizing the I/O performance over the network, it can be advisable to :ref:`annotate the input file access patterns of rules <storage-access-patterns>`.
Snakemake provides lots of tunables for non-local execution, which can all be found under :ref:`all_options` and in the plugin descriptions of the `Snakemake plugin catalog <https://snakemake.github.io/snakemake-plugin-catalog>`__.
In any case, the cluster or cloud specific configuration will entail lots of command line options to be chosen and set, which should be persisted in a :ref:`profile <executing-profiles>`.

---------------------------------
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

.. _executing-profiles:

--------
Profiles
--------

Adapting runs of Snakemake workflows to a particular computing environment can entail many flags and options.
Therefore, since Snakemake 4.1, it is possible to set default options in configuration profile files in YAML format.
Two kinds of profiles are supported:

1. **:ref:`global-profiles`** are used to define default options for a particular system or compute environment, like the default cluster submission command, the default number of jobs to run in parallel or the default amount of memory to reserve for a job.
   They should be applicable to all Snakemake workflows a user runs in that compute environment.
2. A **:ref:`workflow-specific-profile`** profile (introduced in Snakemake 7.29) is used to define default and rule-specific :ref:`snakefiles-resources` specifications for a particular workflow instance. 

.. _profile-files:

Profile YAML files
^^^^^^^^^^^^^^^^^^

The default naming pattern for profile YAML files is ``profile.v9+.yaml``, where the version specifier infix ``v9+.`` is optional.
This naming pattern is required when you refer to profiles by (directory) name or relative path (to directory containing the actual YAML file), or if you want to specify a minimum required version of snakemake via the optional infix (``vX+``).
If you directly reference the actual YAML file by name, you can use an arbitrary name for the profile YAML file.

Alongside the actual profile YAML file, the profile folder can additionally contain auxiliary files.
These can for example be jobscripts or wrappers. 
See https://github.com/snakemake/snakemake-cluster-profiles for examples.

While the different types of profiles should usually contain distinct sets of settings, you can configure any of Snakemake's command line arguments in any of these profiles.
However, if you also provide the same argument in the ``snakemake`` call on the command line, this command line specification will always take precedence.

For example, a :ref:`global-profiles` YAML file with

.. code-block:: yaml

    executor: slurm
    jobs: 110
    default-resources:
      mem_mb: 1024

would set Snakemake to always submit to the SLURM cluster using the respective executor plugin, and to never use more than 110 parallel jobs in total.
It gets interpreted into setting ``--executor slurm --jobs 110 --default-resources mem_mb=1024`` on the command line.

For more complex (nested) options, you can use standard YAML nesting syntax; and for simple switch flags, you can set or unset them with the values ``True`` and ``False``, respectively.
So, for example, this YAML map in a :ref:`workflow-specific-profile`

.. code-block:: yaml

    keep-going: True
    set-threads:
      myrule: 5
    set-resources:
      myrule:
        mem_mb: 500

will be parsed to ``--keep-going --set-threads myrule=5 --set-resources myrule:mem_mb=500``.
Alternatively, you can also specify anything below the top level keys as a string.
So the following would parse to the same command line argument setup:

.. code-block:: yaml

    set-threads: myrule=5
    set-resources: myrule:mem_mb=500

All of these resource specifications can also be made dynamic, by using expressions and certain variables that are available.
For details of the variables you can use, refer to the callable signatures given in the documentation sections on the specification of :ref:`threads <snakefiles-threads>` and :ref:`dynamic resources <snakefiles-dynamic-resources>`.
These enable ``profile.yaml`` entries like:

.. code-block:: yaml

    default-resources:
      mem_mb: max(1.5 * input.size_mb, 100)
    set-threads:
      myrule: max(input.size_mb / 5, 2)
    set-resources:
      myrule:
        mem_mb: attempt * 200


Also, values in profiles can make use of globally available environment variables, for example the ``$USER`` variable.
For example, the following entry would set the default prefix for storing local copies of remote storage files to a user specific directory

.. code-block:: yaml

    local-storage-prefix: /local/work/$USER/snakemake-scratch

Any such environment variables are automatically expanded when evaluating the profile.

Finally, we recommend annotating such profiles with clear comments.
From experience, the most useful mode is usually to include comments right above a setting, including the reasoning behind the chosen value and linkouts to any documentation with further information.
For inspiration, see the examples in the following sections.

.. _global-profiles:

Global profiles
^^^^^^^^^^^^^^^

Global profiles are used to define default options for a particular system or compute environment, applicable to all Snakemake workflows run on a particular system.

.. _defining-global-profiles:

Defining global profiles
""""""""""""""""""""""""

Default options specified in global profiles will include things like the default cluster submission command, the default number of jobs to run in parallel or the default amount of memory to reserve for a job.
We recommend to clearly motivate any configuration choices in comments, for example

.. code-block:: yaml

    # This cluster uses the slurm job submission system, for details on how to
    # configure the respective executor plugin for snakemake, see
    # https://snakemake.github.io/snakemake-plugin-catalog/plugins/executor/slurm.html
    executor: slurm
    # This cluster allows users to run 100 jobs concurrently. As the slurm
    # executor plugin only checks for completed jobs with a lag, we slightly
    # oversubmit jobs to always have jobs available in the queue.
    jobs: 110
    # This cluster allows the use of 200 cores per user. As with the number of
    # jobs, we slightly oversubmit.
    cores: 220
    # If a rule doesn't have the following resources specified, it will default
    # to requesting the resources specified here.
    default-resources:
      mem_mb: 1024

Also, before creating your own global profile from scratch, check whether someone has already created and shared such a profile for your local compute environment at:
https://github.com/snakemake/snakemake-cluster-profiles
This repository can also serve as inspiration, when creating a new profile.
And if you have created such a profile for your local compute cluster, feel free to share it via this repository.
Just make sure to check with your local system administrators, if all the information included is OK to be shared publicly.

.. _using-global-profiles:

Using global profiles
"""""""""""""""""""""

To make use of such a global profile, you can make Snakemake aware of it in one of two ways:

1. You set the environment variable ``$SNAKEMAKE_PROFILE``.
2. You use the command line argument ``--profile``.

In either of these cases, you can reference a profile in multiple ways:

1. With a profile name, which represents a subdirectory of one of the standard locations, and is assumed to contain a file called ``profile.yaml`` (or ``config.yaml`` for backwards compatibility).
2. With a path to such a directory containing a ``profile.yaml`` file, relative to one of the standard locations.
3. With a path to a YAML file with an arbitrary name, relative to one of the standard locations.

The standard locations that snakemake searches are the current working directory, a standard system-wide location and a standard user-specific location.
As the system and user locations are system-dependent, you should always check the ``--profile`` entry of the ``snakemake --help``.
This will list the locations specific to the system that you run the command on.
On Linux, these locations are ``/etc/xdg/snakemake`` (system-wide) and ``$HOME/.config/snakemake`` (user-specific).

As this leads to a plethora of ways to specify profiles, let us provide some examples.

As a first example, assume a system-wide profile with the absolute path ``/etc/xdg/snakemake/system_profile/profile.v9+.yaml`` on a Linux system.
You can set the ``$SNAKEMAKE_PROFILE`` variable or the ``--profile`` argument to any of:

.. code-block:: bash

    system_profile
    system_profile/profile.v9+.yaml
    /etc/xdg/snakemake/system_profile/profile.v9+.yaml

As a second example, assume a user-specific profile with absolute path ``$HOME/.config/snakemake/user_profile/profile.yaml`` on a Linux system.
You can set the ``$SNAKEMAKE_PROFILE`` variable or the ``--profile`` argument to any of:

.. code-block:: bash

    user_profile
    user_profile/profile.yaml
    $HOME/.config/snakemake/user_profile/profile.yaml

As a third example, assume a user-specific profile in a custom location, ``/path/to/user_profile/custom_profile_name.yaml``.
In this case you have to always use the absolute path because the profile is not in a standard location, and you also have to reference the file by name, as it doesn't follow the standard naming convention.

The environment variable can either be set by system administrators, providing a default profile for all users of a system.
Alternatively, users can set this environment variable for themselves.
Usually, this is done by setting the following in your shell startup configuration (for example in your ``~/.bashrc`` for ``bash`` shells):

.. code-block:: bash

    export SNAKEMAKE_PROFILE="~/.config/snakemake/my_profile/profile.v9+.yaml"

If you instead provide the profile via the ``--profile`` command line argument, the ``$SNAKEMAKE_PROFILE`` environment variable will be ignored.

.. _using-multiple-global-profiles:

Using multiple global profiles
""""""""""""""""""""""""""""""

If multiple instances of the ``--profile`` command line argument are given, all the profiles are merged.
While merging, the profile instances specified later take precedence over earlier instances, wherever the same top-level entries occur in multiple profiles.
Take the following example files and their invocation.

``/path/to/system_profile.yaml``

.. code-block:: yaml

    executor: slurm
    cores: 100
    default-resources:
      mem_mb: 1024
      # some resource, where you can only use two at any time
      rate_limiter: 2

``~/.config/snakemake/user_profile/profile.yaml``

.. code-block:: yaml

    cores: 50
    default-resources:
      mem_mb: 8000

When loaded in the order

.. code-block:: bash

    snakemake --profile /path/to/system_profile.yaml --profile user_profile

this will lead to snakemake being run with the configuration

.. code-block:: yaml

    executor: slurm
    cores: 50
    default-resources:
      mem_mb: 8000

Thus, the ``user_profile`` takes precedence over the ``system_profile.yaml``.
As the ``user_profile`` does not specify an ``executor``, ``executor: slurm`` is kept.
As the ``user_profile`` also specifies ``cores``, its entry of ``cores: 50`` overwrites the value from the ``system_profile.yaml``.
And as you can see from the last example, this overwriting of entries happens at the top level of YAML entries.
As the ``user_profile`` specifies the top-level entry ``default-resources``, this whole entry from ``system_profile.yaml`` is discarded and replaced by what is specified in ``user_profile``.
Thus, the ``default-resources: rate_limiter: 2`` entry is lost.

Workflow specific profiles take precedence over global profiles in the same way.
And arguments specified on the command line take precedence over any profile YAML files.

.. _workflow-specific-profile:

Workflow specific profile
^^^^^^^^^^^^^^^^^^^^^^^^^

A workflow specific profile (introduced in Snakemake 7.29) is used to define default and rule-specific :ref:`snakefiles-resources` specifications for a particular workflow. 

.. _defining-workflow-specific-profiles:

Defining workflow specific profiles
"""""""""""""""""""""""""""""""""""

Mostly, a workflow-specific profile is meant to set rule-specific resources and command-line arguments that are specific to the running of that workflow.

.. code-block:: yaml

    # rule-specific threads settings
    set-threads:
      # my_rule is set to threads: 4 in the rule, but we are trying out more
      # parallelism on this instance. If this proves useful, we'll propagate
      # this back to the rule definition.
      my_rule: 8
    # default non-threads resources and any arbitrary resources a workflow
    # defines, can be controlled here
    set-resources:
      my_rule:
        # This particular dataset seems to need more memory, while usually
        # the mem_mb=80000 from the rule is enough. Maybe we can create a
        # dynamic resource allocation based on the following issue: <link>
        mem_mb: 400000
        # This much memory is only available in this dedicated slurm partition.
        slurm_partition: high_memory
    # See docs on defining scatter-gather processes
    set-scatter:
      # For this dataset, the default scatter for this rule of 200 is far too
      # much. We can optimise it by splitting into fewer but bigger chunks of
      # data.
      scatter_rule: 4

But a workflow specific profile can also overwrite ``default-resources`` from :ref:`global-profiles`, for example

.. code-block:: yaml

    # default resources that are assigned if a rule doesn't have mem_mb
    # specified via its `resources:` directive.
    default-resources:
      # While a lower amount of memory might be a good default for other
      # workflows, most rules in this one just need a higher amount.
      mem_mb: 32000

So, as you can already see from the examples, most of the settings in these kinds of profiles should eventually be propagated back into (dynamic) resource settings for every individual rule (via the ``resources:`` directive).
But they are a very good tool for a number of purposes.
For example, to quickly change things on the fly, without having to change anything in the underlying workflow and waiting for another release.
Or, to distribute workflow profiles that optimise the rules' resource usage of a workflow for a particular computing environment.
Especially for the latter, it can be useful to distribute workflow specific profiles along with the workflow itself.
For example, when the workflow has its Snakefile at ``workflow/Snakefile``, a profile tailored to a particular ``xyz_cluster`` could be placed at ``workflow/profiles/xyz_cluster/profile.yaml`` and then used with ``--workflow-profile xyz_cluster``.

.. _using-workflow-specific-profiles:

Using workflow specific profiles
""""""""""""""""""""""""""""""""

To make use of such a workflow specific profile, you can make Snakemake aware of it in one of two ways:

1. You give it a standardised filename (for example ``profile.v9+.yaml``, see the section on :ref:`profile-files`) and save it in the default folder hierarchy relative to the ``Snakefile`` or the current working directory, usually either ``profiles/default/profile.yaml`` or ``workflow/profiles/default/profile.yaml``.
2. You use the command line argument ``--workflow-profile``.

Note that even without specifying ``--workflow-profile``, Snakemake will automatically search for and apply a workflow profile in ``profiles/default/`` (relative to the Snakefile or working directory).
To prevent any workflow profile from being loaded, you can explicitly call ``--workflow-profile none``, as using the command line argument ensures that the ``default/`` location is not searched implicitly (unless you explicitly specify ``--workflow-profile default``).

.. note::

  When using modules, the profile will not be propagated to the main workflow importing that module.
  However, using `snakedeploy deploy-workflow` to deploy a workflow as a module, will also copy any profiles included under the standard location `workflow/profiles` (for more info, see `the snakedeploy documentation for deploying workflows <https://snakedeploy.readthedocs.io/en/stable/workflow_users/workflow_deployment.html>`_).
  Starting from this import, or starting with a new file, users can create a profile for that main workflow.

Any profile you specify on the command line is searched in paths relative to the ``Snakefile`` location and the current working directory.
You have the same options to specify it as for :ref:`using-global-profiles`:

1. With a profile name, which represents a subdirectory of one of the standard locations, and is assumed to contain a file called ``profile.yaml`` (or ``config.yaml`` for backwards compatibility).
2. With a path to such a directory containing a ``profile.yaml`` file, relative to one of the standard locations.
3. With a path to a YAML file with an arbitrary name, relative to one of the standard locations.

For example, if your ``Snakefile`` sits in the recommended location in subfolder ``workflow/``, ``snakemake --workflow-profile my_profile`` will look for:

.. code-block:: bash

    profiles/my_profile/profile.yaml
    workflow/profiles/my_profile/profile.yaml

Note, the examples here omit the optional ``vX+`` minimum version infix.

With the same ``Snakefile`` location, ``snakemake --workflow-profile relative_path/to/my_profile`` will look for:

.. code-block:: bash

    relative_path/to/my_profile/profile.yaml
    profiles/relative_path/to/my_profile/profile.yaml
    workflow/profiles/relative_path/to/my_profile/profile.yaml

And finally, assuming that the specified file exists, ``snakemake --workflow-profile extra_profiles_dir/workflow_profile.yaml`` will short-circuit the lookup and just use the file that is specified.

Whenever a workflow profile is successfully specified, it is parsed after any global profiles.
It takes precedence over them, overriding any pre-existing top-level keys that it also specifies, but keeping any top-level keys that it doesn't contain.

For example, if the ``--profile global_profile`` YAML file sets

.. code-block:: yaml

    cores: 50
    default-resources:
      mem_mb: 8000
      disk_mb: 20000

and the ``--workflow-profile workflow_on_xyz`` sets

.. code-block:: yaml

    default-resources:
      mem_mb: 4000
      extra: something
    keep-going: True

the resulting profile configuration will be

.. code-block:: yaml

    cores: 50
    default-resources:
      mem_mb: 4000
      extra: something
    keep-going: True

Similarly, any specifications in your workflow specific profile will be overwritten by command line arguments of the ``snakemake`` run.
So, if you run ``snakemake --workflow-profile workfklow_on_xyz --default-resources mem_mb 1000 --cores 2``, the resulting configuration will be:

.. code-block:: yaml

    cores: 2
    default-resources:
      mem_mb: 1000
    keep-going: True


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

