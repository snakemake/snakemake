.. _manual-main:

=====================================
Welcome to Snakemake's documentation!
=====================================

.. image:: https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg
    :target: https://bioconda.github.io/recipes/snakemake/README.html

.. image:: https://img.shields.io/pypi/pyversions/snakemake.svg
    :target: https://www.python.org

.. image:: https://img.shields.io/pypi/v/snakemake.svg
    :target: https://pypi.python.org/pypi/snakemake

.. image:: https://quay.io/repository/snakemake/snakemake/status
       :target: https://quay.io/repository/snakemake/snakemake

.. image:: https://img.shields.io/circleci/project/bitbucket/snakemake/snakemake.svg
    :target: https://circleci.com/bb/snakemake/snakemake/tree/master

.. image:: https://img.shields.io/badge/stack-overflow-orange.svg
    :target: http://stackoverflow.com/questions/tagged/snakemake

Snakemake is an MIT-licensed workflow management system that aims to reduce the complexity of creating workflows by providing a fast and comfortable execution environment, together with a clean and modern specification language in python style.
Snakemake workflows are essentially Python scripts extended by declarative code to define **rules**.
Rules describe how to create **output files** from **input files**.


.. _manual-quick_example:

-------------
Quick Example
-------------

.. code-block:: python

    rule targets:
        input:
            "plots/dataset1.pdf",
            "plots/dataset2.pdf"

    rule plot:
        input:
            "raw/{dataset}.csv"
        output:
            "plots/{dataset}.pdf"
        shell:
            "somecommand {input} {output}"


* Similar to GNU Make, you specify targets in terms of a pseudo-rule at the top.
* For each target and intermediate file, you create rules that define how they are created from input files.
* Snakemake determines the rule dependencies by matching file names.
* Input and output files can contain multiple named wildcards.
* Rules can either use shell commands, plain Python code or external Python or R scripts to create output files from input files.
* Snakemake workflows can be executed on workstations and clusters without modification. The job scheduling can be constrained by arbitrary resources like e.g. available CPU cores, memory or GPUs.
* Snakemake can use Amazon S3, Google Storage, Dropbox, FTP and SFTP to access input or output files and further access input files via HTTP and HTTPS.

.. _main-getting-started:

---------------
Getting started
---------------

To get started, consider the :ref:`tutorial <tutorial-welcome>`, the `introductory slides <http://slides.com/johanneskoester/deck-1>`_, and the :ref:`FAQ <project_info-faq>`.

.. _main-support:

-------
Support
-------

* In case of questions, please post on `stack overflow <http://stackoverflow.com/questions/tagged/snakemake>`_.
* To discuss with other Snakemake users, you can use the `mailing list <https://groups.google.com/forum/#!forum/snakemake>`_.
* For bugs and feature requests, please use the `issue tracker <https://bitbucket.org/snakemake/snakemake/issues>`_.
* For contributions, visit Snakemake on `bitbucket <https://bitbucket.org/snakemake/snakemake>`_ and read the :ref:`guidelines <project_info-contributing>`.

--------
Citation
--------

`Köster, Johannes and Rahmann, Sven. "Snakemake - A scalable bioinformatics workflow engine". Bioinformatics 2012. <http://bioinformatics.oxfordjournals.org/content/28/19/2520>`_

See :doc:`Citations <project_info/citations>` for more information.

----------------
Related Projects
----------------

`Snakemake Wrappers Repository <https://snakemake-wrappers.readthedocs.org>`_
    The Snakemake Wrapper Repository is a collection of reusable wrappers that allow to quickly use popular command line tools from Snakemake rules and workflows.

`Snakemake Workflow Repository <https://bitbucket.org/snakemake/snakemake-workflows>`_
    This repository provides a collection of high quality modularized and re-usable rules and workflows.
    The provided code should also serve as a best-practices of how to build production ready workflows with Snakemake.
    Everybody is invited to contribute.

`Bioconda <https://bioconda.github.io/>`_
    Bioconda can be used from Snakemake for creating completely reproducible workflows by pin pointing the used software version and providing binaries.


.. toctree::
   :caption: Installation
   :name: installation
   :hidden:
   :maxdepth: 1

   getting_started/installation
   getting_started/examples


.. toctree::
   :caption: Tutorial
   :name: tutorial
   :hidden:
   :maxdepth: 1

   tutorial/welcome
   tutorial/basics
   tutorial/advanced
   tutorial/additional_features

.. toctree::
  :caption: Executing workflows
  :name: execution
  :hidden:
  :maxdepth: 1

  executable.rst

.. toctree::
    :caption: Defining workflows
    :name: snakefiles
    :hidden:
    :maxdepth: 1

    snakefiles/writing_snakefiles
    snakefiles/rules
    snakefiles/configuration
    snakefiles/modularization
    snakefiles/remote_files
    snakefiles/utils
    snakefiles/deployment


.. toctree::
    :caption: API Reference
    :name: api-reference
    :hidden:
    :maxdepth: 1

    api_reference/snakemake
    api_reference/snakemake_utils


.. toctree::
    :caption: Project Info
    :name: project-info
    :hidden:
    :maxdepth: 1

    project_info/citations
    project_info/more_resources
    project_info/faq
    project_info/contributing
    project_info/authors
    project_info/history
    project_info/license
