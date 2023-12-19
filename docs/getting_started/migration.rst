.. _migration:

====================================
Migration between Snakemake versions
====================================

Snakemake is meant to remain backwards compatible as much as possible.
However, sometimes, very rarely, we remove old almost unused features that have since then
been replaced by new ones (so far, this happened only once, for Snakemake 8).
Sometimes, new features are added that do not require, but make it strongly advisable to adapt workflows (e.g. because the new features provide a better user or recipient experience).

Below are migration hints for particular Snakemake versions.

Migrating to Snakemake 8
------------------------

Workflow definitions
^^^^^^^^^^^^^^^^^^^^

Snakemake 8 removes the support for four syntactical elements, which are all officially deprecated since multiple major releases:

* Support for marking output files as ``dynamic`` has been removed. You should instead use :ref:`checkpoints <snakefiles-checkpoints>`.
* Support for the ``version`` directive has been removed. You should use the :ref:`conda <integrated_package_management>` or :ref:`container <apptainer>` integration instead.
* Support for the ``subworkflow`` directive has been removed. You should use the :ref:`module directive <snakefiles-modules>` instead, which provides the same functionality in a more general way.
* Support for remote providers has been removed. You should use :ref:`storage plugins <storage-support>` instead. 
  Most of the old remote providers have been migrated into the new storage plugins
  (see the `Snakemake plugin catalog <https://snakemake.github.io/snakemake-plugin-catalog>`_.).
  Two former remote providers have been migrated into Snakemake wrappers instead, namely 
  the NCBI and ENA remote providers, which are now replaced by the 
  `entrez/efetch <https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/entrez/efetch.html>`_ and 
  the `ena <https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/ena.html>`_ wrappers.
  As of writing, the Snakemake storage plugin for xrootd (see `here <https://github.com/snakemake/snakemake-storage-plugin-xrootd>`_) does not yet pass the CI tests. Any help would be greatly appreciated.


Command line interface
^^^^^^^^^^^^^^^^^^^^^^

The command line interface of Snakemake 8 has a lot of new options which are best explored using::

    snakemake --help

Morever, some options have been renamed:

* All the execution backends have been moved into plugins. When you used e.g. ``--kubernetes`` and corresponding options before, you should now use ``--executor kubeternes`` and check the `Snakemake plugin catalog <https://snakemake.github.io/snakemake-plugin-catalog/plugins/executor/kubernetes.html>`_ for the new options. The same holds for all other execution backends, see `here <https://snakemake.github.io/snakemake-plugin-catalog/index.html>`_.
* The ``--use-conda`` and ``--use-singularity`` options are deprecated. Instead you should now use ``--software-deployment-method conda`` or ``--software-deployment-method apptainer`` or ``--software-deployment-method conda apptainer`` if you need both.

Profiles
^^^^^^^^

Profiles can now be versioned.
If your profile makes use of settings that are available in version 8 or later, use the filename ``config.v8+.yaml`` for the profile configuration (see :ref:`profiles`).

API
^^^

The Snakemake API has been completely rewritten into a modern `dataclass <https://docs.python.org/3/library/dataclasses.html>`_ based approach.
The traditional central ``snakemake()`` function is gone.
For an example how to use the new API, check out the Snakemake CLI implementation `here <https://github.com/snakemake/snakemake/blob/04ec2c0262b2cb96cbcd7edbbb2596979c1703ae/snakemake/cli.py#L1767>`_.