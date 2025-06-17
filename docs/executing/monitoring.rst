.. _monitoring:

==========
Monitoring
==========

Since Snakemake 9.0, Snakemake supports monitoring workflow execution through logger plugins.
These plugins enable integration with various monitoring systems and services, allowing you to track workflow progress, collect runtime statistics, and analyze execution patterns.

Logger plugins can capture detailed information about workflow execution, including:

* Job start and completion times
* Resource usage statistics
* Error messages and debugging information
* Workflow topology and dependencies


To use a logger plugin, specify it via the ``--logger`` command line option:

.. code-block:: console

    $ snakemake --logger <plugin-name> --cores 4

Multiple logger plugins can be used simultaneously by specifying the ``--logger`` option multiple times.

Available logger plugins and their configuration options can be found in the `Snakemake plugin catalog <https://snakemake.github.io/snakemake-plugin-catalog>`_.
The catalog provides detailed documentation for each plugin, including installation instructions, configuration examples, and usage guidelines.

For details on developing logger plugins, see the documentation in the `Snakemake logger plugin interface repository <https://github.com/snakemake/snakemake-interface-logger-plugins>`_.
