.. _interaction_visualization_reporting_tutorial:

==================================================================
Tutorial: Interaction, Visualization, and Reporting with Snakemake
==================================================================

.. _Snakemake: https://snakemake.github.io
.. _Jupyter: https://jupyter.org
.. _Datavzrd: https://datavzrd.github.io

Via its Jupyter_ :ref:`notebook <snakefiles_notebook-integration>` and :ref:`script <snakefiles-external_scripts>` integration, Snakemake_ offers powerful support for interactive data exploration.
This is particularly useful for the last mile of data analysis, where results obtained from established tools and libraries are summarized and visualized for e.g. presentation in a scientific publication.
Moreover, via its :ref:`reporting <snakefiles-reports>` capabilities, Snakemake can be used to generate reproducible reports that contain the results of a data analysis pipeline, the code that was used to generate them, and the environment in which the pipeline was executed.
The latter integrates well with Datavzrd_, a tool for generating interactive views of tabular data.

In this short tutorial, we bring all these capabilities together in order to demonstrate how Snakemake enables last mile data analysis without loosing reproducibility, adaptability, and transparency.

Setup
-----

Install Snakemake into a new Conda environment as instructed in the :ref:`installation guide <getting_started-installation>`.
In the following, we assume that you conduct all commands within this environment.
If you already have Snakemake installed, make sure that you use the latest stable version.

Next, create an empty working directory at a reasonable place in your file system.
In the following, we assume that you have an open terminal inside of that directory.
If you use visual studio code (recommended), open a new instance inside of that directory and open a terminal via the terminal menu.

Step 1: Obtain example data
---------------------------

First, we need some example data.
We will leverage this step to conduct our first interactive exploration via Snakemake's Jupyter integration.
We create a new directory ``workflow`` and inside of that an empty ``Snakefile``.
In the ``Snakefile``, define the following rule:

.. code-block:: python

    rule download:
        output:
            "data/iris.csv"
        script:
            "scripts/download.py"