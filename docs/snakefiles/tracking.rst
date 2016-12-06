.. _snakefiles-tracking:

===============
Change Tracking
===============

Snakemake tracks changes in files and rules for re-execution after such a change.


.. _snakefiles-version_tracking:

Version Tracking
================

Rules can specify a version that is tracked by Snakemake together with the output files. When the version changes snakemake informs you when using the flag ``--summary`` or ``--list-version-changes``.
The version can be specified by the version directive, which takes a string:

.. code-block:: python

    rule:
        input:   ...
        output:  ...
        version: "1.0"
        shell:   ...

The version can of course also be filled with the output of a shell command, e.g.:

.. code-block:: python

    SOMECOMMAND_VERSION = subprocess.check_output("somecommand --version", shell=True)

    rule:
        version: SOMECOMMAND_VERSION

Alternatively, you might want to use file modification times in case of local scripts:

.. code-block:: python

    SOMECOMMAND_VERSION = str(os.path.getmtime("path/to/somescript"))

    rule:
        version: SOMECOMMAND_VERSION

A re-run can be automated by invoking Snakemake as follows:

.. code-block:: console

    $ snakemake -R `snakemake --list-version-changes`


.. _snakefiles-code_tracking:

Code Tracking
=============

Snakemake tracks the code that was used to create your files.
In combination with ``--summary`` or ``--list-code-changes`` this can be used to see what files may need a re-run because the implementation changed.
Re-run can be automated by invoking Snakemake as follows:

.. code-block:: console

    $ snakemake -R `snakemake --list-code-changes`
