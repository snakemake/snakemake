.. _executors:

======================
Using executor plugins
======================

By default, Snakemake will run all rules locally on the machine Snakemake itself is running on. 
For distributed computing, i.e. on a cluster or in the cloud, Snakemake offers an interface for so-called executor plugins (as of version 8.0.0 and onward). 

To set up Snakemake for distributed execution on your compute environment, first navigate to the `Snakemake plugin catalog <https://snakemake.github.io/snakemake-plugin-catalog/>`_. 
Locate an executor suitable your cluster/cloud environment and follow its respective installation instructions.

.. admonition:: No fitting executor?
   :class: note

   In this case you may check out the `generic cluster executor plugin <https://snakemake.github.io/snakemake-plugin-catalog/plugins/executor/cluster-generic.html>`_ and configure it with your compute environment's submission and status check commands.
   Or - even better - you can use `Snakedeploy <https://snakedeploy.readthedocs.io/en/stable/snakemake_developers/scaffold_snakemake_plugins.html>`_ to set up a base package structure for `creating your own executor plugin <https://snakemake.github.io/snakemake-plugin-catalog/index.html#create-new-plugins>`_. 

