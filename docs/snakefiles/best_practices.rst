.. _snakefiles-best_practices:

==============
Best practices
==============

Care about code quality
-----------------------

Snakemake (>=5.11) comes with a code quality checker (a so called linter), that analyzes your workflow and highlights issues that should be solved in order to follow best practices, achieve maximum readability, and reproducibility.
The linter can be invoked with 

.. code-block:: bash

      snakemake --lint

given that a ``Snakefile`` or ``workflow/Snakefile`` is accessible from your working directory.
It is **highly recommended** to run the linter before publishing any workflow, asking questions on Stack Overflow or filing issues on Github.

Care about code readability
---------------------------

1. There is an automatic formatter for Snakemake workflows, called `Snakefmt <https://github.com/snakemake/snakefmt>`_, which should be applied to any Snakemake workflow before publishing it.
2. Try to keep filenames short (thus easier on the eye), but informative. Avoid mixing of too many special characters (e.g. decide whether to use ``_`` or ``-`` as a separator and do that consistently throughout the workflow).
3. Try to keep Python code like helper functions separate from rules (e.g. in a ``workflow/rules/common.smk`` file). This way, you help non-experts to read the workflow without needing to parse internals that are irrelevant for them. The helper function names should be chosen in a way that makes them sufficiently informative without looking at their content. Also avoid ``lambda`` expressions inside of rules.
4. Use Snakemake's :ref:`semantic helper functions <snakefiles-semantic-helpers>` in order to increase readibility and to avoid the reimplementation of common functionality for aggregation, parameter lookup or path modifications.

Ensure portability
------------------

Annotate all your rules with versioned :ref:`Conda <integrated_package_management>` or :ref:`container <apptainer>` based software environment definitions. This ensures that your workflow utilizes the exactly same isolated software stacks, independently of the underlying system.

Generate interactive reports (for free)
---------------------------------------

Annotate your final results for inclusing into Snakemake's automatic :ref:`interactive reports <snakefiles-reports>` (thereby make sure to use all the features, including categories and labels).
This makes them explorable in a high-level way, while connecting them to the workflow code, parameters, and software stack.

Enable configurability
----------------------

Configuration of a workflow should be handled via :ref:`config files <snakefiles_standard_configuration>` and, if needed, tabular configuration like sample sheets (either via :ref:`Pandas <snakefiles_tabular_configuration>` or :ref:`PEPs <snakefiles-peps>`).
Use such configuration for metadata and experiment information, **not for runtime specific configuration** like threads, resources and output folders.
For those, just rely on Snakemake's CLI arguments like ``--set-threads``, ``--set-resources``, ``--set-default-resources``, and ``--directory``. 
This makes workflows more readable, scalable, and portable.

Avoid duplication of efforts
----------------------------

Make use of `Snakemake wrappers <https://snakemake-wrappers.readthedocs.io>`_ whenever possible. Consider contributing to the wrapper repo whenever you have a rule that reoccurs in at least two of your workflows.

Test your workflow continuously
-------------------------------

When hosting your workflow in a `Github <https://github.com>`_ repository, it is a good idea to add some minimal test data and configure `Github Actions <https://github.com/features/actions>`_ for continuously testing the workflow on each new commit.For this purpose, we provide predefined Github actions for both running tests and linting `here <https://github.com/snakemake/snakemake-github-action>`__, as well as formatting `here <https://github.com/snakemake/snakefmt#github-actions>`__.

Follow the standards
--------------------

1. For publishing and distributing a Snakemake workflow, it is a good idea to stick to a :ref:`standardized folder structure <distribution_and_reproducibility>` that is expected by frequent users of Snakemake. This simplifies the navigation through the codebase and keeps the workflow repository and the working directory clean.
2. The `Snakemake workflow catalog <https://snakemake.github.io/snakemake-workflow-catalog>`_ automatically lists Snakemake workflows hosted on `Github <https://github.com>`_ if they follow certain `rules <https://snakemake.github.io/snakemake-workflow-catalog/?rules=true>`_.
   By complying to these `rules <https://snakemake.github.io/snakemake-workflow-catalog/?rules=true>`_ you can make your workflow more discoverable and even automate its usage documentation (see `"Standardized usage" <https://snakemake.github.io/snakemake-workflow-catalog/?rules=true>`_).