.. _executor_tutorial:

============================
Snakemake Executor Tutorials
============================

.. _cloud executors: https://snakemake.readthedocs.io/en/stable/executing/cloud.html
.. _tutorial: https://snakemake.readthedocs.io/en/stable/tutorial/tutorial.html

This set of tutorials are intended to introduce you to executing `cloud executors`_.
We start with the original Snakemake `tutorial`_ and expand upon it to be run
in different cloud environments. For each run, we show you how to:

 - authenticate with credentials, if required
 - prepare your workspace
 - submit a basic job
 - generate an error and debug

The examples presented in these tutorials come from Bioinformatics.
However, Snakemake is a general-purpose workflow management system for any discipline.
We ensured that no bioinformatics knowledge is needed to understand the tutorial.

.. toctree::
   :maxdepth: 2

   google_lifesciences
   azure_aks
   flux

   
