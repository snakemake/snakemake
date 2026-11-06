.. _snakefiles-debugging_workflows:

===================
Debugging workflows
===================

Debugging workflows can be a challenging task, especially for complex workflows with many rules and dependencies.
Here are some tips and tools that can help you troubleshoot workflows effectively:

Log files
---------

Each Snakemake run will produce a log file in ``.snakemake/log/``, mirroring the information printed to the console when running Snakemake.
This log file can be especially helpful when running Snakemake in a non-interactive environment, e.g. when executing Snakemake as a cluster job or in a container.


Logs of remotely executed jobs
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Depending on the executor you are using, additional log files may be generated for each job.
For the location thereof, please refer to the `documentation of the respective executor plugin <https://snakemake.github.io/snakemake-plugin-catalog/>`_.

Redirecting STDERR of rules
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

When using the ``log`` directive on rules using ``script``, you can add the following snippets to the rules' scripts to redirect their STDERR to the specified log file:

**For Python scripts:**

.. code-block:: python

    import sys
    sys.stderr = open(snakemake.log[0], "w", buffering=1)

Here, the ``buffering=1`` ensures that line buffering is used, so that ``STDERR`` lines are written to the log file whenever a full line is available.
This avoids information not getting printed before throwing an error due to some longer buffering.


**For R scripts:**

.. code-block:: r

    log <- file(snakemake@log[[1]], open="wt")
    sink(log)
    sink(log, type="message")

Also, if you are looking to have proper backtraces even for unexpected errors (errors not properly handled in your code or in a package you load), you can use:

.. code-block:: r

    rlang::global_entrace()

You will need to have the package ``rlang`` installed, but this for example comes with the ``tidyverse``.
For infos on the function, see the `rlang documentation <https://rlang.r-lib.org/reference/global_entrace.html>`_.
Also, this is `not expected to incur a performance reduction <https://github.com/r-lib/rlang/issues/1717#issuecomment-2163180629>`_.



Saving time on DAG building and resolution
------------------------------------------

To quickly debug a particular rule, you can specify the output of that rule as the desired target when running Snakemake.
This will speed up building the DAG, as Snakemake will only resolve the part of the DAG that is necessary to produce the specified output file.
To avoid clashes with command line argument specifications, it is best to provide the desired output file as the first argument right after ``snakemake``:

.. code-block:: bash

      snakemake path/to/desired/output.file <other arguments>


Debugging in an IDE
-------------------

You can step through Snakemake code in any IDE that can attach to a `debugpy <https://github.com/microsoft/debugpy>`_ server, e.g. Visual Studio Code, Zed or Neovim (nvim-dap).
debugpy comes with Snakemake, so you do not need to install it.
You can debug either the main workflow process or the Python scripts of individual jobs.
In both cases, Snakemake waits on a network address (default ``127.0.0.1:5678``) until your IDE's debugger attaches.

Configuring the IDE
^^^^^^^^^^^^^^^^^^^

Point your IDE's debugger at the address that Snakemake waits on.
In Visual Studio Code, add this configuration to ``.vscode/launch.json``:

.. code-block:: json

    {
        "name": "DebugPy: attach to Snakemake",
        "type": "debugpy",
        "request": "attach",
        "connect": {"host": "127.0.0.1", "port": 5678}
    }

You only need to set this once, because the same configuration works for the main workflow process and for scripts.
To attach, select it in the "Run and Debug" tab and start it (keyboard shortcut: ``F5``).

.. warning:: Anyone who can connect to the debugpy port can execute arbitrary code as your user. Keep the host ``127.0.0.1``. To debug on a remote machine, e.g. an HPC login node, forward the port via SSH instead of listening on a public address:

   .. code-block:: console

       ssh -L 5678:127.0.0.1:5678 user@remote-host

Debugging the main workflow process
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Debugging the main process lets you inspect Snakemake's ``workflow`` object, e.g. in :ref:`event handlers <snakefiles-job_lifetime_handlers>`, which helps when developing both workflows and plugins.

1. Add ``breakpoint()`` statements to your Snakefile where execution should halt, e.g. in the ``onerror`` handler:

   .. code-block:: python

       onerror:
           breakpoint()
           print("Workflow failed")

2. Start Snakemake through debugpy, which waits for a debugger to attach before it runs the workflow:

   .. code-block:: console

       python -m debugpy --listen 127.0.0.1:5678 --wait-for-client -m snakemake --cores 1

3. Attach your IDE's debugger.

Snakemake then halts at each ``breakpoint()``, with the local variables of the surrounding Snakefile code available for inspection.

.. note:: The IDE highlights the wrong line in the Snakefile.
   Snakemake translates the Snakefile into Python code before it runs it, and many directives take more lines in that code than in the Snakefile.
   For example, ``onstart:`` becomes two lines, a decorator and a function definition.
   The IDE shows the line number of the translated code, so the highlighted line moves further down the Snakefile with every directive above the breakpoint.
   Execution still halts at the right place, and the variables and call stack that the IDE shows are correct.
   To check where you are, read the function name in the call stack, e.g. ``__onstart`` or ``__onerror``.
   For the same reason, set breakpoints with ``breakpoint()`` statements rather than in the IDE.

This does not debug individual jobs, because Snakemake runs them in separate processes.

.. _snakefiles-debugging_scripts_in_ide:

Debugging Python scripts of individual jobs
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

When executing Snakemake with the ``--debug-scripts`` flag, each job of a rule with a Python ``script:`` directive or a Python-based ``wrapper:`` pauses before the script runs and waits for your IDE's debugger to attach.
Breakpoints that you set in the script file in your IDE are then hit, with the ``snakemake`` object available for inspection.

Requirements and limitations:

* debugpy does not need to be installed in the job's software environment (e.g. a conda environment of the rule).
  However, the job's Python has to be version 3.9 or newer, because debugpy does not support older versions. A job with an older Python fails before the script runs.
* Execution must be local and on a single core (``--cores 1``), so that only one job waits for a debugger at a time.
* Only Python scripts and Python-based wrappers pause. Other job types, e.g. ``run:``, ``shell:``, ``notebook:`` or scripts in other languages, run as usual.
* A job with ``benchmark: repeat(..., n)`` runs once instead of ``n`` times, so that it pauses only once. Its benchmark file includes the time spent waiting for the debugger and at breakpoints.
* Jobs that run in a container (``container:`` with ``--sdm container``) are not supported, because the container does not see debugpy and Snakemake from the host environment.

To debug a script:

1. Start Snakemake, e.g. targeting the output of the rule you want to debug:

   .. code-block:: console

       snakemake --cores 1 --debug-scripts path/to/desired/output.file

2. When a script job starts, Snakemake prints which rule is paused, the script file to set breakpoints in and the address to attach to.
3. Set breakpoints in that script file and attach your IDE's debugger.
4. The script runs once the debugger is attached and stops at your breakpoints. When the script finishes, the debugger disconnects and Snakemake continues. The next script job pauses again, so attach the debugger again for it. Press ``Ctrl-C`` in the terminal to abort a waiting job.

To use another address, set the environment variables ``SNAKEMAKE_DEBUGPY_HOST`` and ``SNAKEMAKE_DEBUGPY_PORT``.
To debug the main workflow process and scripts at the same time, start the main process on another port, e.g. ``--listen 127.0.0.1:5679``, and add a second IDE configuration for that port.


Debugging individual jobs without an IDE
----------------------------------------

**For Python scripts / run blocks:**

When executing Snakemake with the ``--debug`` flag, Snakemake will drop into an interactive `Python debugger <https://docs.python.org/3/library/pdb.html>`_ (PDB) session. By including ``breakpoint()`` statements in your code you can specify where PDB should halt execution, allowing you to explore the current state of the job.

**For R scripts / run blocks:**

You can save the entire current state of a workspace in R for debugging.
Insert this line right before the code that triggers an error:

.. code-block:: r

    save.image(file = "workspace.RData")


Activate the conda environment that the rule uses (you can find this in the rule's logging output, with a statement ``Activating conda environment: <path to environment>``) and start an interactive R session.
In this session, load all the ``library()`` s that you need for the script.
Then you can load the full workspace and interactively explore / debug what's going on:

.. code-block:: r

    load("workspace.RData")


Preserving wrapper scripts
--------------------------

Snakemake produces a series of wrapper scripts for rules using the ``script`` directive (default location ``.snakemake/scripts/``).
Normally, these are deleted after each run. For debugging purposes, you can disable this behavior by running snakemake with the ``--skip-script-cleanup`` flag.
