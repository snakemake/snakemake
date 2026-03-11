.. _provenance:

==========
Provenance
==========

By default, snakemake stores provenance information / metadata in the ``.snakemake/metadata`` directory.
However, for workflows with large numbers of inputs/outputs, this can lead to issues with the underlying filesystem,
especially when using networked filesystems or filesystems with limitations on the number of files in a directory.

To address this, an experimental DB based provenance system can be enabled using

.. code-block:: console

    $ snakemake --persistence-backend db [--persistence-backend-db-url URL]

By default, this will store provenance information in a SQLite database located at ``.snakemake/metadata.db``.
However, users can specify a different database URL using the ``--persistence-backend-db-url`` option, which supports any database backend supported by SQLAlchemy;
note that the database backend must support JSON columns (e.g. PostgreSQL, MySQL, SQLite 3).