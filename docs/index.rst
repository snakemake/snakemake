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

Snakemake is a workflow management system that aims to reduce the complexity of creating workflows by providing a fast and comfortable execution environment, together with a clean and modern specification language in python style.
Snakemake workflows are essentially Python scripts extended by declarative code to define **rules**.
Rules describe how to create **output files** from **input files**.

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


Features
--------

* Similar to GNU Make, you specify targets in terms of a pseudo rule at the top.
* For each target and intermediate file, you create rules that define how they are created from from input files.
* Snakemake determines the rule dependencies by matching file names.
* Input and output files can contain multiple named wildcards.
* Rules can either use shell commands, plain Python code or external Python or R scripts to create output files from input files.
* Snakemake workflows can be executed on workstations and clusters without modification. The job scheduling can be constrained by arbitrary resources like e.g. available CPU cores, memory or GPUs.
* Snakemake can use Amazon S3, Google Storage, Dropbox, FTP, SFTP to access input or output files and further access input files via HTTP and HTTPS.


.. toctree::
   :caption: Installation & Getting Started
   :name: getting-started
   :hidden:
   :maxdepth: 1

   getting_started/installation
   getting_started/getting_started
   getting_started/examples


.. toctree::
   :caption: Tutorial
   :name: tutorial
   :hidden:
   :maxdepth: 1


.. toctree::
    :caption: Snakefile's
    :name: snakefiles
    :hidden:
    :maxdepth: 1

    snakefiles/writing_snakefiles
    snakefiles/grammar
    snakefiles/rules
    snakefiles/remote_files
    snakefiles/targets
    snakefiles/configuration
    snakefiles/cluster_configuration
    snakefiles/lifetime_handlers
    snakefiles/includes
    snakefiles/sub_workflows
    snakefiles/resources
    snakefiles/utils
    snakefiles/reports
    snakefiles/r_scripting
    snakefiles/job_properties
    snakefiles/depend_version
    snakefiles/dynamic_files
    snakefiles/input_functions
    snakefiles/tracking
    snakefiles/wrappers


.. toctree::
    :caption: User Manual
    :name: user-manual
    :hidden:
    :maxdepth: 1

    user_manual/invocation
    user_manual/bash_completion
    user_manual/visualization
    user_manual/all_options


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

    project_info/support
    project_info/faq
    project_info/contributing
    project_info/authors
    project_info/history
    project_info/license


.. toctree::
    :caption: Citing, Citations, etc.
    :hidden:

    cit_art_etc/citations
    cit_art_etc/citing
    cit_art_etc/talks_posters
    cit_art_etc/workflow_repo
    cit_art_etc/external_resources
