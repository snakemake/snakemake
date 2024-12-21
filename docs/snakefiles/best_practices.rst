.. _snakefiles-best_practices:

==============
Best practices
==============

* Snakemake (>=5.11) comes with a code quality checker (a so called linter), that analyzes your workflow and highlights issues that should be solved in order to follow best practices, achieve maximum readability, and reproducibility.
  The linter can be invoked with 

  .. code-block:: bash

      snakemake --lint

  given that a ``Snakefile`` or ``workflow/Snakefile`` is accessible from your working directory.
  It is **highly recommended** to run the linter before publishing any workflow, asking questions on Stack Overflow or filing issues on Github.
* There is an automatic formatter for Snakemake workflows, called `Snakefmt <https://github.com/snakemake/snakefmt>`_, which should be applied to any Snakemake workflow before publishing it.
* When publishing your workflow in a `Github <https://github.com>`_ repository, it is a good idea to add some minimal test data and configure `Github Actions <https://github.com/features/actions>`_ for continuously testing the workflow on each new commit.
  For this purpose, we provide predefined Github actions for both running tests and linting `here <https://github.com/snakemake/snakemake-github-action>`__, as well as formatting `here <https://github.com/snakemake/snakefmt#github-actions>`__.
* For publishing and distributing a Snakemake workflow, it is a good idea to stick to a :ref:`standardized structure <distribution_and_reproducibility>` that is expected by frequent users of Snakemake.
  The `Snakemake workflow catalog <https://snakemake.github.io/snakemake-workflow-catalog>`_ automatically lists Snakemake workflows hosted on `Github <https://github.com>`_ if they follow certain `rules <https://snakemake.github.io/snakemake-workflow-catalog/?rules=true>`_.
  By complying to these `rules <https://snakemake.github.io/snakemake-workflow-catalog/?rules=true>`_ you can make your workflow more discoverable and even automate its usage documentation (see `"Standardized usage" <https://snakemake.github.io/snakemake-workflow-catalog/?rules=true>`_).
* Configuration of a workflow should be handled via :ref:`config files <snakefiles_standard_configuration>` and, if needed, tabular configuration like sample sheets (either via :ref:`Pandas <snakefiles_tabular_configuration>` or :ref:`PEPs <snakefiles-peps>`).
  Use such configuration for metadata and experiment information, **not for runtime specific configuration** like threads, resources and output folders.
  For those, just rely on Snakemake's CLI arguments like ``--set-threads``, ``--set-resources``, ``--set-default-resources``, and ``--directory``. 
  This makes workflows more readable, scalable, and portable.
* Try to keep filenames short (thus easier on the eye), but informative. Avoid mixing of too many special characters (e.g. decide whether to use ``_`` or ``-`` as a separator and do that consistently throughout the workflow).
* Try to keep Python code like helper functions separate from rules (e.g. in a ``workflow/rules/common.smk`` file). This way, you help non-experts to read the workflow without needing to parse internals that are irrelevant for them. The helper function names should be chosen in a way that makes them sufficiently informative without looking at their content. Also avoid ``lambda`` expressions inside of rules.
* Make use of `Snakemake wrappers <https://snakemake-wrappers.readthedocs.io>`_ whenever possible. Consider contributing to the wrapper repo whenever you have a rule that reoccurs in at least two of your workflows.
