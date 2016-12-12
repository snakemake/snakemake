.. user_manual-snakemake_executable:

============================
The ``snakemake`` Executable
============================

This part of the documentation describes the ``snakemake`` executable.  Snakemake
is primarily a command-line tool, so the ``snakemake`` executable is the primary way
to execute, debug, and visualize workflows.

.. user_manual-snakemake_options:

-----------------------------
Useful Command Line Arguments
-----------------------------

If called without parameters, i.e.

.. code-block:: console

    $ snakemake

Snakemake tries to execute the workflow specified in a file called ``Snakefile`` in the same directory (instead, the snakefile can be given via the parameter ``-s``).

By issueing

.. code-block:: console

    $ snakemake -n

a dry-run can be performed.
This is useful to test if the workflow is defined properly and to estimate the amount of needed computation.
Further, the reason for each rule execution can be printed via


.. code-block:: console

    $ snakemake -n -r

Importantly, snakemake can automatically determine which parts of the workflow can be run in parallel.
By specifying the number of available cores, i.e.

.. code-block:: console

    $ snakemake -j 4

one can tell snakemake to use up to 4 cores and solve a binary knapsack problem to optimize the scheduling of jobs.
If the number is omitted (i.e., only ``-j`` is given), the number of used cores is determined as the number of available CPU cores in the machine.

Finally, snakemake can make use of cluster engines that support shell scripts and have access to a common filesystem, (e.g. the Sun Grid Engine).
In this case, snakemake simply needs to be given a submit command that accepts a shell script as first positional argument:

.. code-block:: console

    $ snakemake --cluster qsub -j 32


Here, -j denotes the number of jobs submitted being submitted to the cluster at the same time (here 32).
The cluster command can be decorated with job specific information, e.g.

.. code-block:: console

    $ snakemake --cluster "qsub {threads}"

Thereby, all keywords of a rule are allowed (e.g. params, input, output, threads, priority, ...).
For example, you could encode the expected running time into params:

.. code-block:: python

    rule:
        input:  ...
        output: ...
        params: runtime="4h"
        shell: ...

and forward it to the cluster scheduler:

.. code-block:: console

    $ snakemake --cluster "qsub --runtime {params.runtime}"

If your cluster system supports `DRMAA <http://www.drmaa.org/>`_, Snakemake can make use of that to increase the control over jobs.
E.g. jobs can be cancelled upon pressing Ctrl+C, which is not possible with the generic ``--cluster`` support.
With DRMAA, no ``qsub`` command need to be provided, but system specific arguments can still be given as a string, e.g.

.. code-block:: console

    $ snakemake --drmaa " -q username" -j 32

Note that the string has to contain a leading whitespace.
Else, the arguments will be interpreted as part of the normal Snakemake arguments, and execution will fail.


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


.. _getting_started-all_options:

-----------
All Options
-----------

All command line options can be printed by calling ``snakemake -h``.  

.. code-block:: text

    usage: snakemake [-h] [--snakefile FILE] [--gui [PORT]] [--cores [N]]
                     [--resources [NAME=INT [NAME=INT ...]]]
                     [--config [KEY=VALUE [KEY=VALUE ...]]] [--list]
                     [--list-target-rules] [--directory DIR] [--dryrun]
                     [--printshellcmds] [--dag] [--rulegraph] [--d3dag]
                     [--summary] [--detailed-summary] [--touch] [--keep-going]
                     [--force] [--forceall] [--forcerun TARGET [TARGET ...]]
                     [--prioritize TARGET [TARGET ...]] [--allow-ambiguity]
                     [--cluster CMD] [--drmaa [ARGS]] [--immediate-submit]
                     [--jobscript SCRIPT] [--jobname NAME] [--reason]
                     [--stats FILE] [--nocolor] [--quiet] [--nolock] [--unlock]
                     [--cleanup-metadata [FILE [FILE ...]]] [--rerun-incomplete]
                     [--ignore-incomplete] [--list-version-changes]
                     [--list-code-changes] [--list-input-changes]
                     [--list-params-changes] [--latency-wait SECONDS]
                     [--wait-for-files [FILE [FILE ...]]] [--benchmark-repeats N]
                     [--notemp] [--keep-target-files]
                     [--allowed-rules ALLOWED_RULES [ALLOWED_RULES ...]]
                     [--timestamp] [--greedyness GREEDYNESS] [--print-compilation]
                     [--overwrite-shellcmd OVERWRITE_SHELLCMD] [--debug]
                     [--profile FILE] [--bash-completion] [--version]
                     [target [target ...]]

    Snakemake is a Python based language and execution environment for GNU Make-
    like workflows.

    positional arguments:
      target                Targets to build. May be rules or files.

    optional arguments:
      -h, --help            show this help message and exit
      --snakefile FILE, -s FILE
                            The workflow definition in a snakefile.
      --gui [PORT]          Serve an HTML based user interface to the given port
                            (default: 8000). If possible, a browser window is
                            opened.
      --cores [N], --jobs [N], -j [N]
                            Use at most N cores in parallel (default: 1). If N is
                            omitted, the limit is set to the number of available
                            cores.
      --resources [NAME=INT [NAME=INT ...]], --res [NAME=INT [NAME=INT ...]]
                            Define additional resources that shall constrain the
                            scheduling analogously to threads (see above). A
                            resource is defined as a name and an integer value.
                            E.g. --resources gpu=1. Rules can use resources by
                            defining the resource keyword, e.g. resources: gpu=1.
                            If now two rules require 1 of the resource 'gpu' they
                            won't be run in parallel by the scheduler.
      --config [KEY=VALUE [KEY=VALUE ...]]
                            Set or overwrite values in the workflow config object.
                            The workflow config object is accessible as variable
                            config inside the workflow. Default values can be set
                            by providing a JSON file (see Documentation).
      --list, -l            Show availiable rules in given Snakefile.
      --list-target-rules, --lt
                            Show available target rules in given Snakefile.
      --directory DIR, -d DIR
                            Specify working directory (relative paths in the
                            snakefile will use this as their origin).
      --dryrun, -n          Do not execute anything.
      --printshellcmds, -p  Print out the shell commands that will be executed.
      --dag                 Do not execute anything and print the directed acyclic
                            graph of jobs in the dot language. Recommended use on
                            Unix systems: snakemake --dag | dot | display
      --rulegraph           Do not execute anything and print the dependency graph
                            of rules in the dot language. This will be less
                            crowded than above DAG of jobs, but also show less
                            information. Note that each rule is displayed once,
                            hence the displayed graph will be cyclic if a rule
                            appears in several steps of the workflow. Use this if
                            above option leads to a DAG that is too large.
                            Recommended use on Unix systems: snakemake --ruledag |
                            dot | display
      --d3dag               Print the DAG in D3.js compatible JSON format.
      --summary, -S         Print a summary of all files created by the workflow.
                            The has the following columns: filename, modification
                            time, rule version, status, plan. Thereby rule version
                            contains the versionthe file was created with (see the
                            version keyword of rules), and status denotes whether
                            the file is missing, its input files are newer or if
                            version or implementation of the rule changed since
                            file creation. Finally the last column denotes whether
                            the file will be updated or created during the next
                            workflow execution.
      --detailed-summary, -D
                            Print a summary of all files created by the workflow.
                            The has the following columns: filename, modification
                            time, rule version, input file(s), shell command,
                            status, plan. Thereby rule version contains the
                            versionthe file was created with (see the version
                            keyword of rules), and status denotes whether the file
                            is missing, its input files are newer or if version or
                            implementation of the rule changed since file
                            creation. The input file and shell command columns are
                            selfexplanatory. Finally the last column denotes
                            whether the file will be updated or created during the
                            next workflow execution.
      --touch, -t           Touch output files (mark them up to date without
                            really changing them) instead of running their
                            commands. This is used to pretend that the rules were
                            executed, in order to fool future invocations of
                            snakemake. Fails if a file does not yet exist.
      --keep-going, -k      Go on with independent jobs if a job fails.
      --force, -f           Force the execution of the selected target or the
                            first rule regardless of already created output.
      --forceall, -F        Force the execution of the selected (or the first)
                            rule and all rules it is dependent on regardless of
                            already created output.
      --forcerun TARGET [TARGET ...], -R TARGET [TARGET ...]
                            Force the re-execution or creation of the given rules
                            or files. Use this option if you changed a rule and
                            want to have all its output in your workflow updated.
      --prioritize TARGET [TARGET ...], -P TARGET [TARGET ...]
                            Tell the scheduler to assign creation of given targets
                            (and all their dependencies) highest priority.
                            (EXPERIMENTAL)
      --allow-ambiguity, -a
                            Don't check for ambiguous rules and simply use the
                            first if several can produce the same file. This
                            allows the user to prioritize rules by their order in
                            the snakefile.
      --cluster CMD, -c CMD
                            Execute snakemake rules with the given submit command,
                            e.g. qsub. Snakemake compiles jobs into scripts that
                            are submitted to the cluster with the given command,
                            once all input files for a particular job are present.
                            The submit command can be decorated to make it aware
                            of certain job properties (input, output, params,
                            wildcards, log, threads and dependencies (see the
                            argument below)), e.g.: $ snakemake --cluster 'qsub
                            -pe threaded {threads}'.
      --drmaa [ARGS]        Execute snakemake on a cluster accessed via DRMAA,
                            Snakemake compiles jobs into scripts that are
                            submitted to the cluster with the given command, once
                            all input files for a particular job are present. ARGS
                            can be used to specify options of the underlying
                            cluster system, thereby using the job properties
                            input, output, params, wildcards, log, threads and
                            dependencies, e.g.: --drmaa ' -pe threaded {threads}'.
                            Note that ARGS must be given in quotes and with a
                            leading whitespace.
      --immediate-submit, --is
                            Immediately submit all jobs to the cluster instead of
                            waiting for present input files. This will fail,
                            unless you make the cluster aware of job dependencies,
                            e.g. via: $ snakemake --cluster 'sbatch --dependency
                            {dependencies}. Assuming that your submit script (here
                            sbatch) outputs the generated job id to the first
                            stdout line, {dependencies} will be filled with space
                            separated job ids this job depends on.
      --jobscript SCRIPT, --js SCRIPT
                            Provide a custom job script for submission to the
                            cluster. The default script resides as 'jobscript.sh'
                            in the installation directory.
      --jobname NAME, --jn NAME
                            Provide a custom name for the jobscript that is
                            submitted to the cluster (see --cluster).NAME is
                            "snakejob.{rulename}.{jobid}.sh" per default. The
                            wildcard {jobid} has to be present in the name.
      --reason, -r          Print the reason for each executed rule.
      --stats FILE          Write stats about Snakefile execution in JSON format
                            to the given file.
      --nocolor             Do not use a colored output.
      --quiet, -q           Do not output any progress or rule information.
      --nolock              Do not lock the working directory
      --unlock              Remove a lock on the working directory.
      --cleanup-metadata [FILE [FILE ...]], --cm [FILE [FILE ...]]
                            Cleanup the metadata of given files. That means that
                            snakemake removes any tracked version info, and any
                            marks that files are incomplete.
      --rerun-incomplete, --ri
                            Re-run all jobs the output of which is recognized as
                            incomplete.
      --ignore-incomplete, --ii
                            Ignore any incomplete jobs.
      --list-version-changes, --lv
                            List all output files that have been created with a
                            different version (as determined by the version
                            keyword).
      --list-code-changes, --lc
                            List all output files for which the rule body (run or
                            shell) have changed in the Snakefile.
      --list-input-changes, --li
                            List all output files for which the defined input
                            files have changed in the Snakefile (e.g. new input
                            files were added in the rule definition or files were
                            renamed). For listing input file modification in the
                            filesystem, use --summary.
      --list-params-changes, --lp
                            List all output files for which the defined params
                            have changed in the Snakefile.
      --latency-wait SECONDS, --output-wait SECONDS, -w SECONDS
                            Wait given seconds if an output file of a job is not
                            present after the job finished. This helps if your
                            filesystem suffers from latency (default 5).
      --wait-for-files [FILE [FILE ...]]
                            Wait --latency-wait seconds for these files to be
                            present before executing the workflow. This option is
                            used internally to handle filesystem latency in
                            cluster environments.
      --benchmark-repeats N
                            Repeat a job N times if marked for benchmarking
                            (default 1).
      --notemp, --nt        Ignore temp() declarations. This is useful when
                            running only a part of the workflow, since temp()
                            would lead to deletion of probably needed files by
                            other parts of the workflow.
      --keep-target-files   Do not adjust the paths of given target files relative
                            to the working directory.
      --allowed-rules ALLOWED_RULES [ALLOWED_RULES ...]
                            Only use given rules. If omitted, all rules in
                            Snakefile are used.
      --timestamp, -T       Add a timestamp to all logging output
      --greedyness GREEDYNESS
                            Set the greedyness of scheduling. This value between 0
                            and 1 determines how careful jobs are selected for
                            execution. The default value (1.0) provides the best
                            speed and still acceptable scheduling quality.
      --print-compilation   Print the python representation of the workflow.
      --overwrite-shellcmd OVERWRITE_SHELLCMD
                            Provide a shell command that shall be executed instead
                            of those given in the workflow. This is for debugging
                            purposes only.
      --debug               Print debugging output.
      --profile FILE        Profile Snakemake and write the output to FILE. This
                            requires yappi to be installed.
      --bash-completion     Output code to register bash completion for snakemake.
                            Put the following in your .bashrc (including the
                            accents): `snakemake --bash-completion` or issue it in
                            an open terminal session.
      --version, -v         show program's version number and exit


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
