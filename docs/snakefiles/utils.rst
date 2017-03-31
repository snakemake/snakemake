.. _snakefiles-utils:

=====
Utils
=====

The module ``snakemake.utils`` provides a collection of helper functions for common tasks in Snakemake workflows. Details can be found in :ref:`utils-api`.

.. _snakefiles-reports:

Reports
-------

The ``report`` function provides an easy mechanism to write reports containing your results. A report is written in reStructuredText_ and compiled to HTML. The function allows you to embed your generated tables and plots into the HTML file. By referencing the files from your text, you can easily provide a semantical connection between them. For using this function, you need to have the docutils_ package installed.

.. _reStructuredText: http://docutils.sourceforge.net/rst.html

.. _docutils: https://pypi.python.org/pypi/docutils

.. code-block:: python

    from snakemake.utils import report

    SOMECONSTANT = 42

    rule report:
        input:  F1="someplot.pdf",
                T1="sometable.txt"
        output: html="report.html"
        run:
            report("""
            =======================
            The title of the report
            =======================

            Write your report here, explaining your results. Don't fear to use math
            it will be rendered correctly in any browser using MathJAX,
            e.g. inline :math:`\sum_{{j \in E}} t_j \leq I`,
            or even properly separated:

            .. math::

                |cq_{{0ctrl}}^i - cq_{{nt}}^i| > 0.5

            Include your files using their keyword name and an underscore: F1_, T1_.

            Access your global and local variables like within shell commands, e.g. {SOMECONSTANT}.
            """, output.html, metadata="Johannes Köster (johannes.koester@uni-due.de)", **input)

The optional metadata argument allows to provide arbitrary additional information to the report, e.g. the author name.
The unpacked input files (``**input``) in the report function generates a list of keyword args, that can be referenced inside the document with the mentioned underscore notation. The files will be embedded into the HTML file using `data URLs <http://en.wikipedia.org/wiki/Data_URI_scheme>`_, thus making the report fully portable and not dependent on your local filesystem structure.


.. _snakefiles-r_scripting:

Scripting with R
----------------

The ``R`` function allows you to use R code in your rules. It relies on `rpy2 <https://pypi.python.org/pypi/rpy2>`_:

.. code-block:: python

    from snakemake.utils import R

    SOMECONSTANT = 42

    rule:
        input:  ...
        output: ...
        run:
            R("""
            # write your R code here
            # Access any global or local variables from the Snakefile with the braces notation
            sqrt({SOMECONSTANT});
            # be sure to mask braces used in R control flow by doubling them:
            if(TRUE) {{
                # do something
            }}
            """)

If you compiled your Python installation from source, make sure that Python was build with sqlite support, which is needed for rpy2.
