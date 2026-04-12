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
By default, snakemake automatically detects if the database is located on a network filesystem (like NFS, CephFS, Lustre, or GPFS) and applies specific SQLite3 optimizations (see `sqlite docs <https://sqlite.org/pragma.html>`__):

* ``PRAGMA journal_mode=PERSIST``: Overwrites the journal header with zeros instead of deleting the file.
* ``PRAGMA synchronous=OFF``: Hands data off to the OS immediately without waiting for disk syncs.
* ``PRAGMA temp_store=MEMORY``: Prevents the creation of temporary lock files over the network.
* ``PRAGMA cache_size=-64000``: Allocates up to ~64MB of RAM for the page cache to minimize network reads.

If a non-network filesystem is detected, snakemake uses: ``PRAGMA synchronous=NORMAL`` and ``PRAGMA journal_mode=TRUNCATE``.

By default, we configure a ``PRAGMA busy_timeout={max(10s, latency_wait)}`` (so 10 seconds by default, or the value of ``latency_wait`` if it is higher) to mitigate locking issues with SQLite.

This may or may not be a problem depending on your specific infrastructure.
If you encounter issues with sqlite, consider using a dedicated database server (like PostgreSQL or MySQL) instead.
Since this backend is experimental, finding the optimal setup for your cluster might require some experimentation. Use at your own risk.