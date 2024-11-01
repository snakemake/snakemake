.. _tutorial:

=====================
Tutorial: General use
=====================

.. _Snakemake: https://snakemake.readthedocs.io
.. _GNU Make: https://www.gnu.org/software/make
.. _Python: https://www.python.org
.. _slides: https://slides.com/johanneskoester/snakemake-tutorial
.. _Conda: https://conda.io
.. _Singularity: https://www.sylabs.io

This tutorial introduces the text-based workflow system Snakemake_.
Snakemake follows the `GNU Make`_ paradigm: workflows are defined in terms of rules that define how to create output files from input files.
Dependencies between the rules are determined automatically, creating a DAG (directed acyclic graph) of jobs that can be automatically parallelized.

Snakemake sets itself apart from other text-based workflow systems in the following way.
Hooking into the Python interpreter, Snakemake offers a definition language that is an extension of Python_ with syntax to define rules and workflow specific properties.
This allows Snakemake to combine the flexibility of a plain scripting language with a pythonic workflow definition.
The Python language is known to be concise yet readable and can appear almost like pseudo-code.
The syntactic extensions provided by Snakemake maintain this property for the definition of the workflow.
Further, Snakemake's scheduling algorithm can be constrained by priorities, provided cores and customizable resources and it provides a generic support for distributed computing (e.g., cluster or batch systems).
Hence, a Snakemake workflow scales without modification from single core workstations and multi-core servers to cluster or batch systems.
Finally, Snakemake integrates with the package manager Conda_ and the container engine Singularity_ such that defining the software stack becomes part of the workflow itself.

The examples presented in this tutorial come from Bioinformatics.
However, Snakemake is a general-purpose workflow management system for any discipline.
We ensured that no bioinformatics knowledge is needed to understand the tutorial.

Also have a look at the corresponding slides_.

As follow-up to this tutorial, we recommend to have a look at the :ref:`interaction, visualization and reporting tutorial <interaction_visualization_reporting_tutorial>`, which focuses on Snakemake's ability to cover the last mile of data analysis, i.e., the generation of publication ready reports and interactive visualizations.


.. toctree::
   :maxdepth: 2

   setup
   basics
   advanced
   additional_features
