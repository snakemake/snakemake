.. _Mamba: https://github.com/mamba-org/mamba

.. _storage-support:

===============
Storage support
===============

By default, input and output files or directories defined and used in Snakemake rules are
written to the local filesystem, thereby interpreting relative paths relative to the current working directory.

However, Snakemake also allows to transparently map input and output files to storage providers.
Storage providers are implemented via plugins.
All available storage providers and their documentation are available in the `Snakemake plugin catalog <https://snakemake.github.io/snakemake-plugin-catalog>`_.


Deployment
----------

Storage plugins can be used by deploying them to your local environment via pip::
   
   pip install snakemake-storage-plugin-s3

or Mamba_::

   mamba install -c conda-forge -c bioconda snakemake-storage-plugin-s3

If you choose to register a storage plugin within your workflow (see below), it is advisable to add
the respective plugin package as a dependency of your workflow itself
(see :ref:`global-workflow-dependencies`).

Usage
-----

In general, there are four ways to use a storage provider.

1. Use the ``--default-storage-provider`` command line argument to set a default storage provider.
   This will be used for all input and output files that are not explicitly mapped to a storage provider.
2. Register a storage provider in the workflow and use it only for particular input and output, not all of them.
3. Register multiple entities of the same storage provider with different names/tags. This allows to e.g. use the same protocol to access multiple different remote storages.
4. Let Snakemake automatically find a matching storage provider.

Using the S3 storage plugin, we will provide an example for all of the cases below.
For provider specific options (also for all options of the S3 plugin which are omitted here for brevity) and all available plugins see the `Snakemake plugin catalog <https://snakemake.github.io/snakemake-plugin-catalog>`_.

As default provider
^^^^^^^^^^^^^^^^^^^
If you want all your input and output (which is not explicitly marked to come from 
another storage) to be written to and read from this storage, you can use it as a 
default provider via::

    snakemake --default-storage-provider s3 --default-storage-prefix s3://mybucket/

Custom settings can be passed as well::

    snakemake --default-storage-provider s3 --default-storage-prefix s3://mybucket/ \
        --storage-s3-max-requests-per-second 10


.. _snakefiles-storage-local-files:

Local input/output files
""""""""""""""""""""""""

Despite using a default storage provider, you might have certain files in your workflow
that still come from the local filesystem. Likewise, when importing a module while
specifying a prefix (see :ref:`snakefiles-modules`), you might have some input files
that come from outside the workflow. In either cases, you can use the ``local`` flag::

    rule example:
        input:
            local("resources/example-input.txt")
        output:
            "example-output.txt"
        shell:
            "..."

Here, ``resources/example-input.txt`` will be interpreted as a local file, while
``example-output.txt`` will be written to the default storage provider that you
have specified (with the prefix prepended).

Note that source paths (see :ref:`snakefiles-aux_source_files`) are also not mapped to the default storage provider.
There is no need to additionally mark them as local.

Within the workflow
^^^^^^^^^^^^^^^^^^^

If you want to use this storage plugin only for specific items, you can register it
inside of your workflow::

    # register storage provider (not needed if no custom settings are to be defined here)
    storage:
        provider="s3",
        # optionally add custom settings here if needed
        # alternatively they can be passed via command line arguments
        # starting with --storage-s3-...
        max_requests_per_second=10,

    rule example:
        input:
            storage.s3(
                # define query to the storage backend here
                ...
            ),
        output:
            "example.txt"
        shell:
            "..."

Using multiple entities of the same storage plugin
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In case you have to use this storage plugin multiple times, but with different settings
(e.g. to connect to different storage servers), you can register it multiple times,
each time providing a different tag::

    # register shared settings
    storage:
        provider="s3",
        # optionally add custom settings here if needed
        # alternatively they can be passed via command line arguments
        # starting with --storage-s3-...
        max_requests_per_second=10,

    # register multiple tagged entities
    storage awss3:
        provider="s3",
        endpoint_url="s3.us-east-2.amazonaws.com"

    rule example:
        input:
            storage.awss3(
                # define query to the storage backend here
                ...
            ),
        output:
            "example.txt"
        shell:
            "..."


Automatic inference
^^^^^^^^^^^^^^^^^^^

If the query for a storage plugin is unique given those plugins that you have currently installed,
you can let Snakemake automatically infer the plugin to use::

    rule example:
        input:
            storage("s3://mybucket/example.txt")
        output:
            "example.txt"
        shell:
            "..."

Credentials
^^^^^^^^^^^

Depending on the storage provider, you might have to provide credentials.
Usually, this can be done via environment variables, e.g. for S3::

    export SNAKEMAKE_STORAGE_S3_ACCESS_KEY=...
    export SNAKEMAKE_STORAGE_S3_SECRET_KEY=...
