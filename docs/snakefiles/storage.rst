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

.. _default_storage:

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

Retrieving and keeping remote files locally
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

When the input file is a remote file, the default behaviour is to download the remote
file to the local area and then remove it after the workflow no longer needs it.

This behaviour can be configured from the command line arguments, using::

    snakemake --keep-storage-local-copies --not-retrieve-storage

where ``--keep-storage-local-copies`` directs snakemake to keep the local copies of
remote files that it makes and ``--not-retrieve-storage`` directs snakemake to not download
copies of remote files.

Additionally, this behaviour can be set at the level of the storage directive e.g.::

    storage:
        provider="http",
        retrieve=False,

    storage http_local:
        provider="http",
        keep_local=True,

    rule example_remote:
        input:
            storage.http("http://example.com/example.txt")
        output:
            "example_remote.txt"
        shell:
            "..."

    rule example_local:
        input:
            storage.http_local("http://example.com/example.txt")
        output:
            "example_local.txt"
        shell:
            "..."

Finally, ``retrieve`` and ``keep_local`` can also be set inside the call to the storage
plugin within a rule::

    storage:
        provider="http",
        retrieve=False,

    rule example_remote:
        input:
            storage.http("http://example.com/example.txt", retrieve=False)
        output:
            "example_remote.txt"
        shell:
            "..."

    rule example_local:
        input:
            storage.http("http://example.com/example.txt", keep_local=True)
        output:
            "example_local.txt"
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

You can also pass additional arguments to the storage call within a rule. 
These arguments will be forwarded to the storage provider settings.
For example::

    rule example:
        input:
            storage("s3://mybucket/example.txt", retries=10)
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

.. _storage-access-patterns:

Access pattern annotation
^^^^^^^^^^^^^^^^^^^^^^^^^

Storage providers can automatically optimize the provision of files based on how the files will be accessed by the respective job.
For example, if a file is only read sequentially, the storage provider can avoid downloading it and instead mount or symlink it (depending on the protocol) for ondemand access.
This can be beneficial, in particular if the sequential access involves only a small part of an otherwise large file.
The three access patterns that can be annotated are:

* ``access.sequential``: The file is read sequentially either from start to end or in (potentially disjoint) chunks, but always in order from the start to the end.
* ``access.random``: The file is read in a non-sequential order.
* ``access.multi``: The file is read sequentially, but potentially multiple times in parallel.

Snakemake considers an input file eligible for on-demand provisioning if it is accessed sequentially by one job in parallel.
In all other cases, multi-access, random access, or sequential access by multiple jobs in parallel, the storage provider will download the file to the local filesystem before it is accessed by jobs.
In case no access pattern is annotated (the default), Snakemake will also download the file.

The access patterns can be annotated via flags.
Usually, one would define sequential access as the default pattern (it should usually be the most common pattern in a workflow).
This can be done via the ``inputflags`` directive before defining any rule.
For specific files, the access pattern can be annotated by the respective flags ``access.sequential``, ``access.random``, or ``access.multi``.

.. code-block:: python

    inputflags:
    access.sequential


    rule a:
        input:
            access.random("test1.in")  # expected as local copy (because accessed randomly)
        output:
            "test1.out"
        shell:
            "cmd_b {input} {output}"


    rule b:
        input:
            access.multi("test1.out") # expected as local copy (because accessed multiple times)
        output:
            "test2.{dataset}.out"
        shell:
            "cmd_b {input} {output}"


    rule c:
        input:
            "test2.{dataset}.out" # expected as on-demand provisioning (because accessed sequentially, the default defined above)
        output:
            "test3.{dataset}.out"
        shell:
            "cmd_c {input} {output}"

Note that there is no guarantee that the storage provider makes use of this information, since the possibilities can vary between storage protocols and the development stage of the storage plugin.