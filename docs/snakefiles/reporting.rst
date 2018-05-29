.. _snakefiles-reports:

-------
Reports
-------

From Snakemake 5.1 on, it is possible to automatically generate detailed self-contained HTML reports that encompass runtime statistics, provenance information, workflow topology and results.
For the latter, the Snakefile has to be annotated with additional information.
Each output file that shall be part of the report has to be marked with the ``report`` flag, which optionally points to a caption in `restructured text format <http://docutils.sourceforge.net/rst.html>`_ and allows to define a ``category`` for grouping purposes.
Moreover, a global workflow description can be defined via the ``report`` directive.
Consider the following example:

.. code-block:: python

  report: "report/workflow.rst"


  rule all:
      input:
          ["fig1.svg", "fig2.png"]


  rule c:
      output:
          "test.{i}.out"
      singularity:
          "docker://continuumio/miniconda3:4.4.10"
      conda:
          "envs/test.yaml"
      shell:
          "sleep `shuf -i 1-3 -n 1`; touch {output}"


  rule a:
      input:
          expand("test.{i}.out", i=range(10))
      output:
          report("fig1.svg", caption="report/fig1.rst", category="Step 1")
      shell:
          "sleep `shuf -i 1-3 -n 1`; cp data/fig1.svg {output}"


  rule b:
      input:
          expand("test.{i}.out", i=range(10))
      output:
          report("fig2.png", caption="report/fig2.rst", category="Step 2")
      shell:
          "sleep `shuf -i 1-3 -n 1`; cp data/fig2.png {output}"

As can be seen, we define a global description which is contained in the file ``report/workflow.rst``.
In addition, we mark ``fig1.svg`` and ``fig2.png`` for inclusion into the report, while in both cases specifying a caption text via again referring to a restructured text file.
Note the paths to the ``.rst``-files are interpreted relative to the current Snakefile.
Inside the ``.rst``-files you can use `Jinja2 <http://jinja.pocoo.org>`_ templating to access context information.
In case of the global description, you can access the config dictionary via ``{{ snakemake.config }}``, (e.g., use ``{{ snakemake.config["mykey"] }}`` to access the key ``mykey``).
In case of output files, you can access the same values as available with the :ref:`script directive <snakefiles-external_scripts>`.
Moreover, in every ``.rst`` document, you can link to

* the **Rules** section (with ``Rules_``),
* the **Statistics** section (with ``Statistics_``),
* the **Results** section (with ``Results_``),
* any **category** section (with ``Mycategory_``, while ``Mycategory`` is the name given for the category argument of the report flag). E.g., with above example, you could write ``see `Step 2`_`` in order to link to the section with the results that have been assigned to the category ``Step 2``.
* any **file** marked with the report flag (with ``myfile.txt_``, while ``myfile.txt`` is the basename of the file, without any leading directories). E.g., with above example, you could write ``see fig2.png_`` in order to link to the result in the report document.

For details about the hyperlink mechanism of restructured text see `here <http://docutils.sourceforge.net/docs/user/rst/quickref.html#hyperlink-targets>`_.

To create the report simply run

.. code-block:: bash

    snakemake --report report.html

after your workflow has finished.
All other information contained in the report (e.g. runtime statistics) is automatically collected during creation.
These statistics are obtained from the metadata that is stored in the ``.snakemake`` directory inside your working directory.
The report for above example can be found :download:`here <../../tests/test_report/report.html>`.
The full example source code can be found `here <https://bitbucket.org/snakemake/snakemake/src/master/tests/test_report/>`_.

Note that the report can be restricted to particular jobs and results by specifying targets at the command line, analog to normal Snakemake execution.
For example, with

.. code-block:: bash

    snakemake fig1.svg --report report-short.html

the report contains only ``fig1.svg``.
