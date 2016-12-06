.. _snakefiles-r_scripting:

================
Scripting with R
================


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
