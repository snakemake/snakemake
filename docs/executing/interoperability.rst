
================
Interoperability
================

.. _cwl_export:

----------
CWL export
----------

Snakemake workflows can be exported to `CWL <https://www.commonwl.org/>`_, such that they can be executed in any `CWL-enabled workflow engine <https://www.commonwl.org/#Implementations>`_.
Since, CWL is less powerful for expressing workflows than Snakemake (most importantly Snakemake offers more flexible scatter-gather patterns, since full Python can be used), export works such that every Snakemake job is encoded into a single step in the CWL workflow.
Moreover, every step of that workflow calls Snakemake again to execute the job. The latter enables advanced Snakemake features like scripts, benchmarks and remote files to work inside CWL.
So, when exporting keep in mind that the resulting CWL file can become huge, depending on the number of jobs in your workflow.
To export a Snakemake workflow to CWL, simply run

.. code-block:: console

    $ snakemake --export-cwl workflow.cwl

The resulting workflow will by default use the `Snakemake docker image <https://hub.docker.com/r/snakemake/snakemake>`_ for every step, but this behavior can be overwritten via the CWL execution environment.
Then, the workflow can be executed in the same working directory with, e.g.,

.. code-block:: console

    $ cwltool workflow.cwl

Note that due to limitations in CWL, it seems currently impossible to avoid that all target files (output files of target jobs), are written directly to the workdir, regardless of their relative paths in the Snakefile.

Note that export is impossible in case the workflow contains :ref:`checkpoints <snakefiles-checkpoints>` or output files with absolute paths.