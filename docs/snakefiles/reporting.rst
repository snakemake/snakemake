.. _snakefiles-reports:

-------
Reports
-------

From Snakemake 5.1 on, it is possible to automatically generate detailed self-contained HTML reports that encompass runtime statistics, provenance information, workflow topology and results.
**As an example, the report of the Snakemake rolling paper can be found** `here <https://snakemake.github.io/resources/report.html>`__.

For including results into the report, the Snakefile has to be annotated with additional information.
Each output file that shall be part of the report has to be marked with the ``report`` flag, which optionally points to a caption in `restructured text format <https://docutils.sourceforge.io/docs/user/rst/quickstart.html>`_ and allows to define a ``category`` for grouping purposes.
Moreover, a global workflow description can be defined via the ``report`` directive.
Consider the following example:

.. code-block:: python

  report: "report/workflow.rst"


  rule all:
      input:
          ["fig1.svg", "fig2.png", "testdir"]


  rule c:
      output:
          "test.{i}.out"
      container:
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
          expand("{model}.{i}.out", i=range(10))
      output:
          report("fig2.png", caption="report/fig2.rst", category="Step 2", subcategory="{model}")
      shell:
          "sleep `shuf -i 1-3 -n 1`; cp data/fig2.png {output}"

  rule d:
      output:
          report(
              directory("testdir"), 
              patterns=["{name}.txt"], 
              caption="report/somedata.rst", 
              category="Step 3")
      shell:
          "mkdir {output}; for i in 1 2 3; do echo $i > {output}/$i.txt; done"

As can be seen, we define a global description which is contained in the file ``report/workflow.rst``.
In addition, we mark ``fig1.svg`` and ``fig2.png`` for inclusion into the report, while in both cases specifying a caption text via again referring to a restructured text file.
Note the paths to the ``.rst``-files are interpreted relative to the current Snakefile.

Inside the ``.rst``-files you can use `Jinja2 <https://jinja.palletsprojects.com>`_ templating to access context information.
In case of the global description, you can access the config dictionary via ``{{ snakemake.config }}``, (e.g., use ``{{ snakemake.config["mykey"] }}`` to access the key ``mykey``).
In case of output files, you can access the same values as available with the :ref:`script directive <snakefiles-external_scripts>` (e.g., ``snakemake.wildcards``).

When marking files for inclusion in the report, a ``category`` and a ``subcategory`` can be given, allowing to group results in of the report.
For both, wildcards (like ``{model}`` see rule b in the example), are automatically replaced with the respective values from the corresponding job.

The last rule ``d`` creates a directory with several files, here mimicking the case that it is impossible to specify exactly which files will be created while writing the workflow (e.g. it might depend on the data).
Nevertheless, it is still possible to include those files one by one into the report by defining inclusion patterns (here ``patterns=["{name}.txt"]``) along with the report flag.
When creating the report, Snakemake will scan the directory for files matching the given patterns and include all of them in the report.
Wildcards in those patterns are made available in the jinja-templated caption document along with the rules wildcards in the ``snakemake.wildcards`` object.

If the output of a rule is a directory with an HTML file hierarchy, it is also possible to specify an entry-point HTML file for inclusion into the report, instead of the ``patterns`` approach from above.
This works as follows:

.. code-block:: python

    rule generate_html_hierarchy:
        output:
            report(directory("test"), caption="report/caption.rst", htmlindex="test.html")
        shell:
            """
            # mimic writing of an HTML hierarchy
            mkdir test
            cp template.html test/test.html
            mkdir test/js
            echo \"alert('test')\" > test/js/test.js
            """

Defining file labels
~~~~~~~~~~~~~~~~~~~~~

In addition to category, and subcategory, it is possible to define a dictionary of labels for each report item.
By that, the actual filename will be hidden in the report and instead a table with the label keys as columns and the values in the respective row for the file will be displayed.
This can lead to less technical reports that abstract away the fact that the results of the analysis are actually files.
Consider the following modification of rule ``b`` from above:

.. code-block:: python

    rule b:
      input:
          expand("{model}.{i}.out", i=range(10))
      output:
          report(
              "fig2.png", 
              caption="report/fig2.rst", 
              category="Step 2", 
              subcategory="{model}",
              labels={
                  "model": "{model}",
                  "figure": "some plot"
              }
          )
      shell:
          "sleep `shuf -i 1-3 -n 1`; cp data/fig2.png {output}"


Determining category, subcategory, and labels dynamically via functions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Similar to e.g. with input file and parameter definition (see :ref:`snakefiles-input_functions`), ``category`` and a ``subcategory`` and ``labels`` can be specified by pointing to a function that takes ``wildcards`` as the first argument (and optionally in addition ``input``, ``output``, ``params`` in any order).
The function is expected to return a string or number (int, float, numpy types), or, in case of labels, a dict with strings as keys and strings or numbers as values.


Linking between items
~~~~~~~~~~~~~~~~~~~~~

In every ``.rst`` document, you can link to

* the **Workflow** panel (with ``Rules_``),
* the **Statistics** panel (with ``Statistics_``),
* any **category** panel (with ``Mycategory_``, while ``Mycategory`` is the name given for the category argument of the report flag). E.g., with above example, you could write ``see `Step 2`_`` in order to link to the section with the results that have been assigned to the category ``Step 2``.
* any **file** marked with the report flag (with ``myfile.txt_``, while ``myfile.txt`` is the basename of the file, without any leading directories). E.g., with above example, you could write ``see fig2.png_`` in order to link to the result in the report document.

For details about the hyperlink mechanism of restructured text see `here <https://docutils.sourceforge.io/docs/user/rst/quickref.html#hyperlink-targets>`__.

Rendering reports
~~~~~~~~~~~~~~~~~

To create the report simply run

.. code-block:: bash

    snakemake --report report.html

after your workflow has finished.
All other information contained in the report (e.g. runtime statistics) is automatically collected during creation.
These statistics are obtained from the metadata that is stored in the ``.snakemake`` directory inside your working directory.


You can define an institute specific stylesheet with:

.. code-block:: bash

    snakemake --report report.html --report-stylesheet custom-stylesheet.css

In particular, this allows you to e.g. set a logo at the top (by using CSS to inject a background for the placeholder ``<div id="brand">``, or overwrite colors.
For an example custom stylesheet defining the logo, see :download:`here <../../tests/test_report/custom-stylesheet.css>`.
The report for above example can be found :download:`here <../../tests/test_report/expected-results/report.html>` (with a custom branding for the University of Duisburg-Essen).
The full example source code can be found `here <https://github.com/snakemake/snakemake/tree/main/tests/test_report/>`__.

Note that the report can be restricted to particular jobs and results by specifying targets at the command line, analog to normal Snakemake execution.
For example, with

.. code-block:: bash

    snakemake fig1.svg --report report-short.html

the report contains only ``fig1.svg``.
