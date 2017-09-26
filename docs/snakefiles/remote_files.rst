.. _snakefiles-remote_files:

============
Remote files
============

In versions ``snakemake>=3.5``.

The ``Snakefile`` supports a wrapper function, ``remote()``, indicating a file is on a remote storage provider (this is similar to ``temp()`` or ``protected()``). In order to use all types of remote files, the Python packages ``boto``, ``moto``, ``filechunkio``, ``pysftp``, ``dropbox``, ``requests``, ``ftputil``, ``XRootD``, and ``biopython`` must be installed.

During rule execution, a remote file (or object) specified is downloaded to the local ``cwd``, within a sub-directory bearing the same name as the remote provider. This sub-directory naming lets you have multiple remote origins with reduced likelihood of name collisions, and allows Snakemake to easily translate remote objects to local file paths. You can think of each local remote sub-directory as a local mirror of the remote system. The ``remote()`` wrapper is mutually-exclusive with the ``temp()`` and ``protected()`` wrappers.

Snakemake includes the following remote providers, supported by the corresponding classes:

* Amazon Simple Storage Service (AWS S3): ``snakemake.remote.S3``
* Google Cloud Storage (GS): ``snakemake.remote.GS``
* File transfer over SSH (SFTP): ``snakemake.remote.SFTP``
* Read-only web (HTTP[S]): ``snakemake.remote.HTTP``
* File transfer protocol (FTP): ``snakemake.remote.FTP``
* Dropbox: ``snakemake.remote.dropbox``
* XRootD: ``snakemake.remote.XRootD``
* GenBank / NCBI Entrez: ``snakemake.remote.NCBI``
* WebDAV: ``snakemake.remote.webdav``
* GFAL: ``snakemake.remote.gfal``


Amazon Simple Storage Service (S3)
==================================

This section describes usage of the S3 RemoteProvider, and also provides an intro to remote files and their usage.

It is important to note that you must have credentials (``access_key_id`` and ``secret_access_key``) which permit read/write access. If a file only serves as input to a Snakemake rule, read access is sufficient. You may specify credentials as environment variables or in the file ``=/.aws/credentials``, prefixed with ``AWS_*``, as with a standard `boto config <http://boto.readthedocs.org/en/latest/boto_config_tut.html>`_. Credentials may also be explicitly listed in the ``Snakefile``, as shown below:

For the Amazon S3 and Google Cloud Storage providers, the sub-directory used must be the bucket name.

Using remote files is easy (AWS S3 shown):

.. code-block:: python

    from snakemake.remote.S3 import RemoteProvider as S3RemoteProvider
    S3 = S3RemoteProvider(access_key_id="MYACCESSKEY", secret_access_key="MYSECRET")

    rule all:
        input:
            S3.remote("bucket-name/file.txt")

Expand still works as expected, just wrap the expansion:


.. code-block:: python

    from snakemake.remote.S3 import RemoteProvider as S3RemoteProvider
    S3 = S3RemoteProvider()

    rule all:
        input:
            S3.remote(expand("bucket-name/{letter}-2.txt", letter=["A", "B", "C"]))

It is possible to use S3-compatible storage by specifying a different endpoint address as the `host` kwarg in the provider, as the kwargs used in instantiating the provider are passed in to `boto <https://boto.readthedocs.org/en/latest/ref/s3.html#boto.s3.connection.S3Connection>`_:

.. code-block:: python

    from snakemake.remote.S3 import RemoteProvider as S3RemoteProvider
    S3 = S3RemoteProvider(access_key_id="MYACCESSKEY", secret_access_key="MYSECRET", host="mystorage.example.com")

    rule all:
        input:
            S3.remote("bucket-name/file.txt")

Only remote files needed to satisfy the DAG build are downloaded for the workflow. By default, remote files are downloaded prior to rule execution and are removed locally as soon as no rules depend on them. Remote files can be explicitly kept by setting the ``keep_local=True`` keyword argument:

.. code-block:: python

    from snakemake.remote.S3 import RemoteProvider as S3RemoteProvider
    S3 = S3RemoteProvider(access_key_id="MYACCESSKEY", secret_access_key="MYSECRET")

    rule all:
        input: S3.remote('bucket-name/prefix{split_id}.txt', keep_local=True)

If you wish to have a rule to simply download a file to a local copy, you can do so by declaring the same file path locally as is used by the remote file:

.. code-block:: python

    from snakemake.remote.S3 import RemoteProvider as S3RemoteProvider
    S3 = S3RemoteProvider(access_key_id="MYACCESSKEY", secret_access_key="MYSECRET")

    rule all:
        input:
            S3.remote("bucket-name/out.txt")
        output:
            "bucket-name/out.txt"
        run:
            shell("cp {output[0]} ./")

In some cases the rule can use the data directly on the remote provider, in these cases ``stay_on_remote=True`` can be set to avoid downloading/uploading data unnecessarily. Additionally, if the backend supports it, any potentially corrupt output files will be removed from the remote. The default for ``stay_on_remote`` and ``keep_local`` can be configured by setting these properties on the remote provider object:

.. code-block:: python

    from snakemake.remote.S3 import RemoteProvider as S3RemoteProvider
    S3 = S3RemoteProvider(access_key_id="MYACCESSKEY", secret_access_key="MYSECRET", keep_local=True, stay_on_remote=True)

The remote provider also supports a new ``glob_wildcards()`` (see :ref:`glob-wildcards`) which acts the same as the local version of ``glob_wildcards()``, but for remote files:

.. code-block:: python

    from snakemake.remote.S3 import RemoteProvider as S3RemoteProvider
    S3 = S3RemoteProvider(access_key_id="MYACCESSKEY", secret_access_key="MYSECRET")
    S3.glob_wildcards("bucket-name/{file_prefix}.txt")

    # (result looks just like as if the local glob_wildcards() function were used on a locally with a folder called "bucket-name")

Google Cloud Storage (GS)
=========================

Usage of the GS provider is the same as the S3 provider.
For authentication, one simply needs to login via the ``gcloud`` tool before
executing Snakemake, i.e.:

.. code-block:: console

    $ gcloud auth application-default login

In the Snakefile, no additional authentication information has to be provided:

.. code-block:: python

    from snakemake.remote.GS import RemoteProvider as GSRemoteProvider
    GS = GSRemoteProvider()

    rule all:
        input:
            GS.remote("bucket-name/file.txt")

File transfer over SSH (SFTP)
=============================

Snakemake can use files on remove servers accessible via SFTP (i.e. most \*nix servers).
It uses `pysftp <https://pysftp.readthedocs.org/en/release_0.2.8/pysftp.html#pysftp.Connection>`_ for the underlying support of SFTP, so the same connection options exist.
Assuming you have SSH keys already set up for the server you are using in the ``Snakefile``, usage is simple:


.. code-block:: python

    from snakemake.remote.SFTP import RemoteProvider
    SFTP = RemoteProvider()

    rule all:
        input:
            SFTP.remote("example.com/path/to/file.bam")

The remote file addresses used must be specified with the host (domain or IP address) and the absolute path to the file on the remote server. A port may be specified if the SSH daemon on the server is listening on a port other than 22, in either the ``RemoteProvider`` or in each instance of ``remote()``:

.. code-block:: python

    from snakemake.remote.SFTP import RemoteProvider
    SFTP = RemoteProvider(port=4040)

    rule all:
        input:
            SFTP.remote("example.com/path/to/file.bam")

.. code-block:: python


    from snakemake.remote.SFTP import RemoteProvider
    SFTP = RemoteProvider()

    rule all:
        input:
            SFTP.remote("example.com:4040/path/to/file.bam")

The standard keyword arguments used by `pysftp <https://pysftp.readthedocs.org/en/release_0.2.8/pysftp.html#pysftp.Connection>`_ may be provided to the RemoteProvider to specify credentials (either password or private key):

.. code-block:: python

    from snakemake.remote.SFTP import RemoteProvider
    SFTP = RemoteProvider(username="myusername", private_key="/Users/myusername/.ssh/particular_id_rsa")

    rule all:
        input:
            SFTP.remote("example.com/path/to/file.bam")

.. code-block:: python

    from snakemake.remote.SFTP import RemoteProvider
    SFTP = RemoteProvider(username="myusername", password="mypassword")

    rule all:
        input:
            SFTP.remote("example.com/path/to/file.bam")

If you share credentials between servers but connect to one on a different port, the alternate port may be specified in the ``remote()`` wrapper:

.. code-block:: python

    from snakemake.remote.SFTP import RemoteProvider
    SFTP = RemoteProvider(username="myusername", password="mypassword")

    rule all:
        input:
            SFTP.remote("some-example-server-1.com/path/to/file.bam"),
            SFTP.remote("some-example-server-2.com:2222/path/to/file.bam")

There is a ``glob_wildcards()`` function:

.. code-block:: python

    from snakemake.remote.SFTP import RemoteProvider
    SFTP = RemoteProvider()
    SFTP.glob_wildcards("example.com/path/to/{sample}.bam")

Read-only web (HTTP[s])
=======================

Snakemake can access web resources via a read-only HTTP(S) provider.
This provider can be helpful for including public web data in a workflow.

Web addresses must be specified without protocol, so if your URI looks like this:

.. code-block:: text

    http://server3.example.com/path/to/myfile.tar.gz

The URI used in the ``Snakefile`` must look like this:

.. code-block:: text

    server3.example.com/path/to/myfile.tar.gz

It is straightforward to use the HTTP provider to download a file to the `cwd`:

.. code-block:: python

    import os
    from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider

    HTTP = HTTPRemoteProvider()

    rule all:
        input:
            HTTP.remote("www.example.com/path/to/document.pdf", keep_local=True)
        run:
            outputName = os.path.basename(input[0])
            shell("mv {input} {outputName}")

To connect on a different port, specify the port as part of the URI string:

.. code-block:: python

    from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
    HTTP = HTTPRemoteProvider()

    rule all:
        input:
            HTTP.remote("www.example.com:8080/path/to/document.pdf", keep_local=True)

By default, the HTTP provider always uses HTTPS (TLS). If you need to connect to a resource with regular HTTP (no TLS), you must explicitly include ``insecure`` as a ``kwarg`` to ``remote()``:

.. code-block:: python

    from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
    HTTP = HTTPRemoteProvider()

    rule all:
        input:
            HTTP.remote("www.example.com/path/to/document.pdf", insecure=True, keep_local=True)

If the URI used includes characters not permitted in a local file path, you may include them as part of the ``additional_request_string`` in the ``kwargs`` for ``remote()``. This may also be useful for including additional parameters you don not want to be part of the local filename (since the URI string becomes the local file name).

.. code-block:: python

    from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
    HTTP = HTTPRemoteProvider()

    rule all:
        input:
            HTTP.remote("example.com/query.php", additional_request_string="?range=2;3")

If the file requires authentication, you can specify a username and password for HTTP Basic Auth with the Remote Provider, or with each instance of `remote()`.
For different types of authentication, you can pass in a Python ```requests.auth`` object (see `here <http://docs.python-requests.org/en/latest/api/#authentication>`_) the `auth` ``kwarg``.

.. code-block:: python

    from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
    HTTP = HTTPRemoteProvider(username="myusername", password="mypassword")

    rule all:
        input:
            HTTP.remote("example.com/interactive.php", keep_local=True)

.. code-block:: python

    from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
    HTTP = HTTPRemoteProvider()

    rule all:
        input:
            HTTP.remote("example.com/interactive.php", username="myusername", password="mypassword", keep_local=True)

.. code-block:: python

    from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
    HTTP = HTTPRemoteProvider()

    rule all:
        input:
            HTTP.remote("example.com/interactive.php", auth=requests.auth.HTTPDigestAuth("myusername", "mypassword"), keep_local=True)

Since remote servers do not present directory contents uniformly, ``glob_wildcards()`` is __not__ supported by the HTTP provider.

File Transfer Protocol (FTP)
============================

Snakemake can work with files stored on regular FTP.
Currently supported are authenticated FTP and anonymous FTP, excluding FTP via TLS.

Usage is similar to the SFTP provider, however the paths specified are relative to the FTP home directory (since this is typically a chroot):

.. code-block:: python

    from snakemake.remote.FTP import RemoteProvider as FTPRemoteProvider

    FTP = FTPRemoteProvider(username="myusername", password="mypassword")

    rule all:
        input:
            FTP.remote("example.com/rel/path/to/file.tar.gz")

The port may be specified in either the provider, or in each instance of `remote()`:

.. code-block:: python

    from snakemake.remote.FTP import RemoteProvider as FTPRemoteProvider

    FTP = FTPRemoteProvider(username="myusername", password="mypassword", port=2121)

    rule all:
        input:
            FTP.remote("example.com/rel/path/to/file.tar.gz")

.. code-block:: python

    from snakemake.remote.FTP import RemoteProvider as FTPRemoteProvider

    FTP = FTPRemoteProvider(username="myusername", password="mypassword")

    rule all:
        input:
            FTP.remote("example.com:2121/rel/path/to/file.tar.gz")

Anonymous download of FTP resources is possible:

.. code-block:: python

    from snakemake.remote.FTP import RemoteProvider as FTPRemoteProvider
    FTP = FTPRemoteProvider()

    rule all:
        input:
            # only keeping the file so we can move it out to the cwd
            FTP.remote("example.com/rel/path/to/file.tar.gz", keep_local=True)
        run:
            shell("mv {input} ./")

``glob_wildcards()``:

.. code-block:: python

    from snakemake.remote.FTP import RemoteProvider as FTPRemoteProvider
    FTP = FTPRemoteProvider(username="myusername", password="mypassword")

    print(FTP.glob_wildcards("example.com/somedir/{file}.txt"))

Setting `immediate_close=True` allows the use of a large number of remote FTP input files in a job where the endpoint server limits the number of concurrent connections. When `immediate_close=True`, Snakemake will terminate FTP connections after each remote file action (`exists()`, `size()`, `download()`, `mtime()`, etc.). This is in contrast to the default behavior which caches FTP details and leaves the connection open across actions to improve performance (closing the connection upon job termination).  :

.. code-block:: python

    from snakemake.remote.FTP import RemoteProvider as FTPRemoteProvider
    FTP = FTPRemoteProvider()

    rule all:
        input:
            # only keep the file so we can move it out to the cwd
            # This server limits the number of concurrent connections so we need to have Snakemake close each after each FTP action.
            FTP.remote(expand("ftp.example.com/rel/path/to/{file}", file=large_list), keep_local=True, immediate_close=True)
        run:
            shell("mv {input} ./")

``glob_wildcards()``:

.. code-block:: python

    from snakemake.remote.FTP import RemoteProvider as FTPRemoteProvider
    FTP = FTPRemoteProvider(username="myusername", password="mypassword")

    print(FTP.glob_wildcards("example.com/somedir/{file}.txt"))

Dropbox
=======

The Dropbox remote provider allows you to upload and download from your `Dropbox <https://www.dropbox.com>`_ account without having the client installed on your machine. In order to use the provider you  first need to register an "app" on the `Dropbox developer website <https://www.dropbox.com/developers/apps/create>`_, with access to the Full Dropbox. After registering, generate an OAuth2 access token. You will need the token to use the Snakemake Dropbox remote provider.

Using the Dropbox provider is straightforward:

.. code-block:: python

    from snakemake.remote.dropbox import RemoteProvider as DropboxRemoteProvider
    DBox = DropboxRemoteProvider(oauth2_access_token="mytoken")

    rule all:
        input:
            DBox.remote("path/to/input.txt")

``glob_wildcards()`` is supported:

.. code-block:: python

    from snakemake.remote.dropbox import RemoteProvider as DropboxRemoteProvider
    DBox = DropboxRemoteProvider(oauth2_access_token="mytoken")

    DBox.glob_wildcards("path/to/{title}.txt")

Note that Dropbox paths are case-insensitive.

XRootD
=======

Snakemake can be used with `XRootD <http://xrootd.org/>`_ backed storage provided the python bindings are installed.
This is typically most useful when combined with the ``stay_on_remote`` flag to minimise local storage requirements.
This flag can be overridden on a file by file basis as described in the S3 remote. Additionally ``glob_wildcards()`` is supported:

.. code-block:: python

    from snakemake.remote.XRootD import RemoteProvider as XRootDRemoteProvider

    XRootD = XRootDRemoteProvider(stay_on_remote=True)
    file_numbers = XRootD.glob_wildcards("root://eospublic.cern.ch//eos/opendata/lhcb/MasterclassDatasets/D0lifetime/2014/mclasseventv2_D0_{n}.root")

    rule all:
        input:
            XRootD.remote(expand("local_data/mclasseventv2_D0_{n}.root", n=file_numbers))

    rule make_data:
        input:
            XRootD.remote("root://eospublic.cern.ch//eos/opendata/lhcb/MasterclassDatasets/D0lifetime/2014/mclasseventv2_D0_{n}.root")
        output:
            'local_data/mclasseventv2_D0_{n}.root'
        shell:
            'xrdcp {input[0]} {output[0]}'

GenBank / NCBI Entrez
=====================

Snakemake can directly source input files from `GenBank <https://www.ncbi.nlm.nih.gov/genbank/>`_ and other `NCBI Entrez databases <https://www.ncbi.nlm.nih.gov/books/NBK25497/table/chapter2.T._entrez_unique_identifiers_ui/?report=objectonly>`_ if the Biopython library is installed.

.. code-block:: python

    from snakemake.remote.NCBI import RemoteProvider as NCBIRemoteProvider
    NCBI = NCBIRemoteProvider(email="someone@example.com") # email required by NCBI to prevent abuse

    rule all:
        input:
            "size.txt"

    rule download_and_count:
        input:
            NCBI.remote("KY785484.1.fasta", db="nuccore")
        output:
            "size.txt"
        run:
            shell("wc -c {input} > {output}")

The output format and source database of a record retrieved from GenBank is inferred from the file extension specified. For example, ``NCBI.RemoteProvider().remote("KY785484.1.fasta", db="nuccore")`` will download a FASTA file while ``NCBI.RemoteProvider().remote("KY785484.1.gb", db="nuccore")`` will download a GenBank-format file. If the options are ambiguous, Snakemake will raise an exception and inform the user of possible format choices. To see available formats, consult the `Entrez EFetch documentation <https://www.ncbi.nlm.nih.gov/books/NBK25499/table/chapter4.T._valid_values_of__retmode_and/?report=objectonly>`_. To view the valid file extensions for these formats, access ``NCBI.RemoteProvider()._gb.valid_extensions``, or instantiate an ``NCBI.NCBIHelper()`` and access ``NCBI.NCBIHelper().valid_extensions`` (this is a property).

When used in conjunction with ``NCBI.RemoteProvider().search()``, Snakemake and ``NCBI.RemoteProvider().remote()`` can be used to find accessions by query and download them:

.. code-block:: python

    from snakemake.remote.NCBI import RemoteProvider as NCBIRemoteProvider
    NCBI = NCBIRemoteProvider(email="someone@example.com") # email required by NCBI to prevent abuse

    # get accessions for the first 3 results in a search for full-length Zika virus genomes
    # the query parameter accepts standard GenBank search syntax
    query = '"Zika virus"[Organism] AND (("9000"[SLEN] : "20000"[SLEN]) AND ("2017/03/20"[PDAT] : "2017/03/24"[PDAT])) '
    accessions = NCBI.search(query, retmax=3)

    # give the accessions a file extension to help the RemoteProvider determine the
    # proper output type.
    input_files = expand("{acc}.fasta", acc=accessions)

    rule all:
        input:
            "sizes.txt"

    rule download_and_count:
        input:
            # Since *.fasta files could come from several different databases, specify the database here.
            # if the input files are ambiguous, the provider will alert the user with possible options
            # standard options like "seq_start" are supported
            NCBI.remote(input_files, db="nuccore", seq_start=5000)

        output:
            "sizes.txt"
        run:
            shell("wc -c {input} > sizes.txt")

Normally, all accessions for a query are returned from ``NCBI.RemoteProvider.search()``. To truncate the results, specify ``retmax=<desired_number>``. Standard Entrez `fetch query options <https://www.ncbi.nlm.nih.gov/books/NBK25499/#chapter4.EFetch>`_ are supported as kwargs, and may be passed in to ``NCBI.RemoteProvider.remote()`` and ``NCBI.RemoteProvider.search()``.

WebDAV
======

WebDAV support is currently ``experimental`` and available in Snakemake 4.0 and later.

Snakemake supports reading and writing WebDAV remote files. The protocol defaults to ``https://``, but insecure connections
can be used by specifying ``protocol=="http://"``. Similarly, the port defaults to 443, and can be overridden by specifying ``port=##`` or by including the port as part of the file address.

.. code-block:: python

    from snakemake.remote import webdav

    webdav = webdav.RemoteProvider(username="test", password="test", protocol="http://")

    rule a:
        input:
            webdav.remote("example.com:8888/path/to/input_file.csv"),
        shell:
            # do something


GFAL
=====

GFAL support is available in Snakemake 4.1 and later.

Snakemake supports reading and writing remote files via the `GFAL <https://dmc.web.cern.ch/projects/gfal-2/home>`_ command line client (gfal-* commands).
By this, it supports various grid storage protocols like `GridFTP <https://en.wikipedia.org/wiki/GridFTP>`_.
In general, if you are able to use the `gfal-*` commands directly, Snakemake support for GFAL will work as well.

.. code-block:: python

    from snakemake.remote import gfal

    gfal = gfal.RemoteProvider(retry=5)

    rule a:
        input:
            gfal.remote("gridftp.grid.sara.nl:2811/path/to/infile.txt")
        output:
            gfal.remote("gridftp.grid.sara.nl:2811/path/to/outfile.txt")
        shell:
            # do something

Authentication has to be setup in the system, e.g. via certificates in the ``.globus`` directory.
Usually, this is already the case and no action has to be taken.
The keyword argument to the remote provider allows to set the number of retries (10 per default) in case of failed commands (the GRID is usually relatively unreliable).
The latter may be unsupported depending on the system configuration.

Note that GFAL support used together with the flags ``--no-shared-fs`` and ``--default-remote-provider`` enables you
to transparently use Snakemake in a grid computing environment without a shared network filesystem.


Remote cross-provider transfers
===============================

It is possible to use Snakemake to transfer files between remote providers (using the local machine as an intermediary), as long as the sub-directory (bucket) names differ:

.. code-block:: python

    from snakemake.remote.GS import RemoteProvider as GSRemoteProvider
    from snakemake.remote.S3 import RemoteProvider as S3RemoteProvider

    GS = GSRemoteProvider(access_key_id="MYACCESSKEYID", secret_access_key="MYSECRETACCESSKEY")
    S3 = S3RemoteProvider(access_key_id="MYACCESSKEYID", secret_access_key="MYSECRETACCESSKEY")

    fileList, = S3.glob_wildcards("source-bucket/{file}.bam")
    rule all:
        input:
            GS.remote( expand("destination-bucket/{file}.bam", file=fileList) )
    rule transfer_S3_to_GS:
        input:
            S3.remote( expand("source-bucket/{file}.bam", file=fileList) )
        output:
            GS.remote( expand("destination-bucket/{file}.bam", file=fileList) )
        run:
            shell("cp {input} {output}")
