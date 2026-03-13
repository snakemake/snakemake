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
For the location thereof, please refer to the documentation of the respective executor.

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

1. Add breakpoints to your Snakefile by adding ``breakpoint()`` statements wherever you want to halt execution to inspect the state of the workflow.
2. Start the debugger in Visual Studio Code by navigating to the "Run and Debug" tab and selecting the "DebugPy: debug Snakemake workflow" configuration (keyboard shortcut: ``F5``). This will open an interactive debugging console, allowing you to step through the workflow execution and inspect variables..


.. note:: This will **not** work for debugging the execution of individual jobs, regardless of whether they are executed locally or remotely. For this, you can use Snakemake's ``--debug`` flag, see below.


Further reading: `Python debugging in VS Code <https://code.visualstudio.com/docs/python/debugging>`_


Debugging of individual jobs
----------------------------------------

**For Python scripts / run blocks:**

When executing Snakemake with the ``--debug`` flag, Snakemake will drop into an interactive `Python debugger <https://docs.python.org/3/library/pdb.html>`_ (PDB) session. By including ``breakpoint()`` statements in your code you can specify where PDB should halt execution, allowing you to explore the current state of the job.

**For R scripts / run blocks:**

You can save the entire current state of a workspace in R, for debugging you can insert this line right before some code triggers an error:

.. code-block:: r

    save.image(file = "my_dump.RData")


In an interactive R session, first load all the ``library()`` s that you need for the script. Then you can load the full workspace and interactively explore / debug what's going on:

.. code-block:: r

    load("my_dump.RData")


Preserving wrapper scripts
--------------------------

Snakemake produces a series of wrapper scripts for rules using the ``script`` directive (default location ``.snakemake/scripts/``). Normally, these are deleted after each run. For debugging purposes, you can disable this behavior by running snakemake with the ``--skip-script-cleanup`` flag.
