.. _manual-main:

=====================================
Welcome to Snakemake's documentation!
=====================================

.. image:: https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat-square
    :target: https://bioconda.github.io/recipes/snakemake/README.html

.. image:: https://img.shields.io/pypi/pyversions/snakemake.svg?style=flat-square
    :target: https://www.python.org

.. image:: https://img.shields.io/pypi/v/snakemake.svg?style=flat-square
    :target: https://pypi.python.org/pypi/snakemake

.. image:: https://img.shields.io/docker/pulls/johanneskoester/snakemake.svg?style=flat-square
    :target: https://hub.docker.com/r/johanneskoester/snakemake/

.. image:: https://img.shields.io/badge/stack-overflow-orange.svg?style=flat-square
    :target: http://stackoverflow.com/questions/tagged/snakemake

Snakemake is a MIT-licensed workflow management system that aims to reduce the complexity of creating workflows by providing a fast and comfortable execution environment, together with a clean and modern specification language in python style.
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


* Similar to GNU Make, you specify targets in terms of a pseudo rule at the top.
* For each target and intermediate file, you create rules that define how they are created from from input files.
* Snakemake determines the rule dependencies by matching file names.
* Input and output files can contain multiple named wildcards.
* Rules can either use shell commands, plain Python code or external Python or R scripts to create output files from input files.
* Snakemake workflows can be executed on workstations and clusters without modification. The job scheduling can be constrained by arbitrary resources like e.g. available CPU cores, memory or GPUs.
* Snakemake can use Amazon S3, Google Storage, Dropbox, FTP, SFTP to access input or output files and further access input files via HTTP and HTTPS.

--------
Citation
--------

`Köster, Johannes and Rahmann, Sven. "Snakemake - A scalable bioinformatics workflow engine". Bioinformatics 2012. <http://bioinformatics.oxfordjournals.org/content/28/19/2520>`_

See :doc:`Citations <project_info/citations>` for more information.

----------------
Related Projects
----------------

`Snakemake Workflow Repository <https://bitbucket.org/snakemake/snakemake-workflows>`_
    This repository provides a collection of high quality modularized and re-usable rules and workflows.
    The provided code should also serve as a best-practices of how to build production ready workflows with Snakemake.
    Everybody is invited to contribute.

`Snakemake Wrappers Repository <https://bitbucket.org/snakemake/snakemake-wrappers>`_
    The Snakemake Wrapper Repository is a collection of reusable wrappers that allow to quickly use popular command line tools from Snakemake rules and workflows.

`Bioconda <https://bioconda.github.io/>`_
    Bioconda can be used from Snakemake for creating completely reproducible workflows by pin pointing the used software version and providing binaries.

.. _main-support:

-------
Support
-------

* In case of questions, please post on `stack overflow <http://stackoverflow.com/questions/tagged/snakemake>`_.
* To discuss with other Snakemake users, you can use the `mailing list <https://groups.google.com/forum/#!forum/snakemake>`_.
* For bugs and feature requests, please use the `issue tracker <https://bitbucket.org/snakemake/snakemake/issues>`_.


.. toctree::
   :caption: Installation & Getting Started
   :name: getting-started
   :hidden:
   :maxdepth: 1

   getting_started/installation
   getting_started/examples
   getting_started/executable


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
    :caption: Snakefile's
    :name: snakefiles
    :hidden:
    :maxdepth: 1

    snakefiles/writing_snakefiles
    snakefiles/rules
    snakefiles/configuration
    snakefiles/modularization
    snakefiles/reports
    snakefiles/remote_files
    snakefiles/r_scripting


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
