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


Interactive debugging
---------------------

Debugging the main workflow process using Visual Studio Code
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The `Visual Studio Code <https://code.visualstudio.com/>`_ editor comes with a powerful debugger for Python.
With a little tinkering, you can also use it to debug Snakemake's main process, which can be very useful for developing both workflows and plugins for Snakemake.

.. note:: This is neither officially supported by Snakemake nor by Visual Studio Code, and may frequently crash!

To set this up, perform the following steps:

1. Install Visual Studio Code and the `Python debugger extension <https://marketplace.visualstudio.com/items?itemName=ms-python.debugpy>`_
2. Open your workflow in Visual Studio Code and create or open the file ``.vscode/launch.json``.
3. Add the following configuration to the ``launch.json`` file:

.. code-block:: json

    {
        "configurations": [

            // ...

            {
                "name": "DebugPy: debug Snakemake workflow",
                "request": "launch",
                "type": "debugpy",
                "cwd": "${workspaceFolder}",
                "args": [
                    "snakemake",
                    "--debug",
                    "--cores",
                    "1",
                    "--nolock",
                    "--forceall",
                    "--executor",
                    "local",
                ],
                "program": "-m",
                "python": "${command:python.interpreterPath}",
                "console": "internalConsole",
                "redirectOutput": true,
                "internalConsoleOptions": "openOnSessionStart",
                // Don't set justMyCode to 'true' - otherwise breakpoints will be skipped.
                // Technically they do not occur within your code, but within Snakemake's workflow.py
                "justMyCode": false,
            },
        ]
    }

To now debug your workflow:

1. Add breakpoints to your Snakefile by adding ``breakpoint()`` statements wherever you want to halt execution to inspect the state of the workflow - for example, in the :ref:`onerror handler <snakefiles-job_lifetime_handlers>`:

.. code-block:: python

    onerror: # Is executed if the pipeline fails
        breakpoint()

2. Start the debugger in Visual Studio Code by navigating to the "Run and Debug" tab and selecting the "DebugPy: debug Snakemake workflow" configuration (keyboard shortcut: ``F5``).
   This will open an interactive debugging console, allowing you to step through the workflow execution and inspect variables, e.g. Snakemake's ``workflow`` object.


.. note:: This will **not** work for debugging the execution of individual jobs, regardless of whether they are executed locally or remotely. For this, you can use Snakemake's ``--debug`` or ``--debug-scripts`` flags, see below.


Further reading: `Python debugging in VS Code <https://code.visualstudio.com/docs/python/debugging>`_


Debugging of individual jobs
----------------------------------------

**For Python scripts / run blocks in the terminal:**

When executing Snakemake with the ``--debug`` flag, Snakemake will drop into an interactive `Python debugger <https://docs.python.org/3/library/pdb.html>`_ (PDB) session. By including ``breakpoint()`` statements in your code you can specify where PDB should halt execution, allowing you to explore the current state of the job.

.. _snakefiles-debugging_scripts_in_ide:

**For Python scripts in an IDE:**

When executing Snakemake with the ``--debug-scripts`` flag, each job of a rule with a Python ``script:`` directive or a Python-based ``wrapper:`` pauses before the script runs and waits for an IDE debugger to attach via `debugpy <https://github.com/microsoft/debugpy>`_.
Breakpoints that you set in the script file in your IDE are then hit, with the ``snakemake`` object available for inspection.
This works with any IDE that can attach to a debugpy server, e.g. Visual Studio Code, Zed or Neovim (nvim-dap).

Requirements and limitations:

* ``debugpy`` must be installed in the environment that runs Snakemake. It does not need to be installed in the job's software environment (e.g. a conda environment of the rule).
* Execution must be local and on a single core (``--cores 1``), so that only one job waits for a debugger at a time.
* Only Python scripts and Python-based wrappers pause. Other job types, e.g. ``run:``, ``shell:``, ``notebook:`` or scripts in other languages, run as usual.

To debug a script:

1. Start Snakemake, e.g. targeting the output of the rule you want to debug:

   .. code-block:: console

       snakemake --cores 1 --debug-scripts path/to/desired/output.file

2. When a script job starts, Snakemake prints which rule is paused, the script file to set breakpoints in and the address to attach to (default ``127.0.0.1:5678``).
3. Set breakpoints in that script file and attach your IDE's debugger. In Visual Studio Code, use a configuration like this in ``.vscode/launch.json``:

   .. code-block:: json

       {
           "name": "DebugPy: attach to Snakemake script",
           "type": "debugpy",
           "request": "attach",
           "connect": {"host": "127.0.0.1", "port": 5678}
       }
   You only need to set this once; this configuration will be reused for all future debug sessions in your working directory.

4. The script runs once the debugger is attached and stops at your breakpoints. When the script finishes, the debugger disconnects and Snakemake continues. The next script job pauses again, so attach the debugger again for it. Press ``Ctrl-C`` in the terminal to abort a waiting job.

To use another address, set the environment variables ``SNAKEMAKE_DEBUGPY_HOST`` and ``SNAKEMAKE_DEBUGPY_PORT``.

.. warning:: Anyone who can connect to the debugpy port can execute arbitrary code as your user. Keep the default host ``127.0.0.1``. To debug on a remote machine, e.g. an HPC login node, forward the port via SSH instead of listening on a public address:

   .. code-block:: console

       ssh -L 5678:127.0.0.1:5678 user@remote-host

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
