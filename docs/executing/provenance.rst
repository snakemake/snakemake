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

If using an SQLite database on networked filesystems:
Note that this can sometimes experience lock contention and latency issues during highly parallel cluster execution.
This may or may not be a problem depending on your specific infrastructure.
If you encounter database locking errors, consider using a dedicated database server (like PostgreSQL or MySQL) instead.
Since this backend is experimental, finding the optimal setup for your cluster might require some experimentation. Use at your own risk.

By default, we configure a ``PRAGMA busy_timeout={max(10s, latency_wait)}`` (so 10 seconds by default, or the value of ``latency_wait`` if it is higher) to mitigate locking issues with SQLite.