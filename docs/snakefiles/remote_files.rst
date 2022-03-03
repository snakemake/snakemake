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
* Microsoft Azure Blob Storage: ``snakemake.remote.AzBlob``
* File transfer over SSH (SFTP): ``snakemake.remote.SFTP``
* Read-only web (HTTP[S]): ``snakemake.remote.HTTP``
* File transfer protocol (FTP): ``snakemake.remote.FTP``
* Dropbox: ``snakemake.remote.dropbox``
* XRootD: ``snakemake.remote.XRootD``
* GenBank / NCBI Entrez: ``snakemake.remote.NCBI``
* WebDAV: ``snakemake.remote.webdav``
* GFAL: ``snakemake.remote.gfal``
* GridFTP: ``snakemake.remote.gridftp``
* iRODS: ``snakemake.remote.iRODS``
* EGA: ``snakemake.remote.EGA``
* AUTO: an automated remote selector

Amazon Simple Storage Service (S3)
==================================

This section describes usage of the S3 RemoteProvider, and also provides an intro to remote files and their usage.

It is important to note that you must have credentials (``access_key_id`` and ``secret_access_key``) which permit read/write access. If a file only serves as input to a Snakemake rule, read access is sufficient. You may specify credentials as environment variables or in the file ``~/.aws/credentials``, prefixed with ``AWS_*``, as with a standard `boto config <https://boto.readthedocs.org/en/latest/boto_config_tut.html>`_. Credentials may also be explicitly listed in the ``Snakefile``, as shown below:

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

If the AWS CLI is installed it is possible to configure your keys globally. This removes the necessity of hardcoding the keys in the Snakefile. The interactive AWS credentials setup can be done using the following command:

.. code-block:: python

    aws configure

S3 then can be used without the keys.

.. code-block:: python

    from snakemake.remote.S3 import RemoteProvider as S3RemoteProvider
    S3 = S3RemoteProvider()

Finally, it is also possible to overwrite the S3 host via adding a ``host`` argument (taking a URL string) to ``S3RemoteProvider``.

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


Microsoft Azure Blob Storage
=============================

Usage of the Azure Blob Storage provider is similar to the S3 provider. For
authentication, an account name and shared access signature (SAS) or key can be used. If these
variables are not passed directly to AzureRemoteProvider (see
[BlobServiceClient
class](https://docs.microsoft.com/en-us/python/api/azure-storage-blob/azure.storage.blob.blobserviceclient?view=azure-python)
for naming), they will be read from environment variables, named
`AZ_BLOB_ACCOUNT_URL` and `AZ_BLOB_CREDENTIAL`. `AZ_BLOB_ACCOUNT_URL` takes the form
`https://<accountname>.blob.core.windows.net` and may also contain a SAS. If
a SAS is not part of the URL, `AZ_BLOB_CREDENTIAL` has to be set to the SAS or alternatively to
the storage account key.

When using AzBlob as default remote provider you will almost always want to
pass these environment variables on to the remote execution environment (e.g.
Kubernetes) with `--envvars`, e.g
`--envvars AZ_BLOB_ACCOUNT_URL AZ_BLOB_CREDENTIAL`.

.. code-block:: python

    from snakemake.remote.AzBlob import RemoteProvider as AzureRemoteProvider
    AS = AzureRemoteProvider()# assumes env vars AZ_BLOB_ACCOUNT_URL and possibly AZ_BLOB_CREDENTIAL are set

    rule a:
        input:
            AS.remote("path/to/file.txt")




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

If you need to create the output directories in the remote server, you can specify ``mkdir_remote=True``  in the ``RemoteProvider`` constructor.

.. code-block:: python

   from snakemake.remote.SFTP import RemoteProvider
   SFTP = RemoteProvider(mkdir_remote=True)

   rule all:
       input:
           "/home/foo/bar.txt"
       output:
           SFTP.remote('example.com/home/foo/create/dir/bar.txt')
       shell:
           "cp {input} {output}"

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

    https://server3.example.com/path/to/myfile.tar.gz

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
For different types of authentication, you can pass in a Python ```requests.auth`` object (see `here <https://requests.readthedocs.io/en/master/api/#authentication>`_) the `auth` ``kwarg``.

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

Snakemake can be used with `XRootD <https://xrootd.slac.stanford.edu/>`_ backed storage provided the python bindings are installed.
This is typically most useful when combined with the ``stay_on_remote`` flag to minimise local storage requirements.
This flag can be overridden on a file by file basis as described in the S3 remote. Additionally ``glob_wildcards()`` is supported:

.. code-block:: python

    from snakemake.remote.XRootD import RemoteProvider as XRootDRemoteProvider

    XRootD = XRootDRemoteProvider(stay_on_remote=True)
    file_numbers = XRootD.glob_wildcards("root://eospublic.cern.ch//eos/opendata/lhcb/MasterclassDatasets/D0lifetime/2014/mclasseventv2_D0_{n}.root").n

    rule all:
        input:
            expand("local_data/mclasseventv2_D0_{n}.root", n=file_numbers)

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
====

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
For an example see the `surfsara-grid configuration profile <https://github.com/Snakemake-Profiles/surfsara-grid>`_.

GridFTP
=======

GridFTP support is available in Snakemake 4.3.0 and later.

As a more specialized alternative to the GFAL remote provider, Snakemake provides a `GridFTP <https://en.wikipedia.org/wiki/GridFTP>`_ remote provider.
This provider only supports the GridFTP protocol. Internally, it uses the `globus-url-copy <http://toolkit.globus.org/toolkit/docs/latest-stable/gridftp/user/#globus-url-copy>`_ command for downloads and uploads, while all other tasks are delegated to the GFAL remote provider.

.. code-block:: python

    from snakemake.remote import gridftp

    gridftp = gridftp.RemoteProvider(retry=5)

    rule a:
        input:
            gridftp.remote("gridftp.grid.sara.nl:2811/path/to/infile.txt")
        output:
            gridftp.remote("gridftp.grid.sara.nl:2811/path/to/outfile.txt")
        shell:
            # do something

Authentication has to be setup in the system, e.g. via certificates in the ``.globus`` directory.
Usually, this is already the case and no action has to be taken.
The keyword argument to the remote provider allows to set the number of retries (10 per default) in case of failed commands (the GRID is usually relatively unreliable).
The latter may be unsupported depending on the system configuration.

Note that GridFTP support used together with the flags ``--no-shared-fs`` and ``--default-remote-provider`` enables you
to transparently use Snakemake in a grid computing environment without a shared network filesystem.
For an example see the `surfsara-grid configuration profile <https://github.com/Snakemake-Profiles/surfsara-grid>`_.


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


iRODS
=====

You can access an iRODS server to retrieve data from and upload data to it.
If your iRODS server is not set to a certain timezone, it is using UTC. It is
advised to shift the modification time provided by iRODS (``modify_time``)
then to your timezone by providing the ``timezone`` parameter such that
timestamps coming from iRODS are converted to the correct time.

iRODS actually does not save the timestamp from your original file but creates
its own timestamp of the upload time. When iRODS downloads the file for
processing, it does not take the timestamp from the remote file. Instead,
the file will have the timestamp when it was downloaded. To get around this,
we create a metadata entry to store the original file stamp from your system
and alter the timestamp of the downloaded file accordingly. While uploading,
the metadata entries ``atime``, ``ctime`` and ``mtime`` are added. When this
entry does not exist (because this module didn't upload the file), we fall back
to the timestamp provided by iRODS with the above mentioned strategy.

To access the iRODS server you need to have an iRODS environment configuration
file available and in this file the authentication needs to be configured.
The iRODS configuration file can be created by following the `official
instructions
<https://docs.irods.org/master/system_overview/configuration/#irodsirods_environmentjson>`_).

The default location for the configuration file is
``~/.irods/irods_environment.json``.  The ``RemoteProvider()`` class accepts
the parameter ``irods_env_file`` where an alternative path to the
``irods_environment.json`` file can be specified.  Another way is to export the
environment variable ``IRODS_ENVIRONMENT_FILE`` in your shell to specify the
location.

There are several ways to configure the authentication against the iRODS
server, depending on what your iRODS server offers. If you are using the
authentication via password, the default location of the authentication file is
``~/.irods/.irodsA``. Usually this file is generated with the ``iinit`` command
from the ``iCommands`` program suite. Inside the ``irods_environment.json``
file, the parameter ``"irods_authentication_file"`` can be set to specifiy an
alternative location for the ``.irodsA`` file. Another possibility to change
the location is to export the environment variable
``IRODS_AUTHENTICATION_FILE``.

The ``glob_wildcards()`` function is supported.

.. code-block:: python

    from snakemake.remote.iRODS import RemoteProvider

    irods = RemoteProvider(irods_env_file='setup-data/irods_environment.json',
                           timezone="Europe/Berlin") # all parameters are optional

    # please note the comma after the variable name!
    # access: irods.remote(expand('home/rods/{f}), f=files))
    files, = irods.glob_wildcards('home/rods/{files})

    rule all:
        input:
            irods.remote('home/rods/testfile.out'),

    rule gen:
        input:
            irods.remote('home/rods/testfile.in')
        output:
            irods.remote('home/rods/testfile.out')
        shell:
            r"""
            touch {output}
            """

An example for the iRODS configuration file (``irods_environment.json``):

.. code-block:: json

    {
        "irods_host": "localhost",
        "irods_port": 1247,
        "irods_user_name": "rods",
        "irods_zone_name": "tempZone",
        "irods_authentication_file": "setup-data/.irodsA"
    }


Please note that the ``zone`` folder is not included in the path as it will be
taken from the configuration file. The path also must not start with a ``/``.

By default, temporarily stored local files are removed. You can specify anyway
the parameter ``overwrite`` to tell iRODS to overwrite existing files that are
downloaded, because iRODS complains if a local file already exists when a
download attempt is issued (uploading is not a problem, though).

In the Snakemake source directory in ``snakemake/tests/test_remote_irods`` you
can find a working example.


EGA
===

The European Genome-phenome Archive (EGA) is a service for permanent archiving
and sharing of all types of personally identifiable genetic and phenotypic data
resulting from biomedical research projects.

From version 5.2 on, Snakemake provides experimental support to use EGA as a remote provider, such that
EGA hosted files can be transparently used as input.
For this to work, you need to define your username and password as environment
variables ``EGA_USERNAME`` and ``EGA_PASSWORD``.

Files in a dataset are addressed via the pattern ``ega/<dataset_id>/<filename>``.
Note that the filename should not include the ``.cip`` ending that is sometimes displayed in EGA listings:

.. code-block:: python

    import snakemake.remote.EGA as EGA

    ega = EGA.RemoteProvider()


    rule a:
        input:
            ega.remote("ega/EGAD00001002142/COLO_829_EPleasance_TGENPipe.bam.bai")
        output:
            "data/COLO_829BL_BCGSC_IlluminaPipe.bam.bai"
        shell:
            "cp {input} {output}"

Upon download, Snakemake will automatically decrypt the file and check the MD5 hash.

Zenodo
======

`Zenodo <https://zenodo.org>`_ is a catch-all open data and software repository. 
Snakemake allows file upload and download from Zenodo. 
To access your Zenodo files you need to set up Zenodo account and create a personal access token with at least write scope.
Personal access token must be supplied as ``access_token`` argument.
You need to supply deposition id as ``deposition`` to upload or download files from your deposition.
If no deposition id is supplied, Snakemake creates a new deposition for upload.
Zenodo UI and REST API responses were designed with having in mind uploads of a total of 20-30 files.
Avoid creating uploads with too many files, and instead group and zip them to make it easier their distribution to end-users.

.. code-block:: python
    from snakemake.remote.zenodo import RemoteProvider
    import os

    # let Snakemake assert the presence of the required environment variable
    envvars:
        "MYZENODO_PAT"

    access_token=os.environ["MYZENODO_PAT"]
    zenodo = RemoteProvider(deposition="your deposition id", access_token=access_token)

    rule upload:
        input:
            "output/results.csv"
        output:
            zenodo.remote("results.csv")
        shell:
            "cp {input} {output}"


It is possible to use `Zenodo sandbox environment <https://sandbox.zenodo.org>`_ for testing by setting ``sandbox=True`` argument.
Using sandbox environment requires setting up sandbox account with its personal access token.

Auto
====

A wrapper which automatically selects an appropriate remote provider based on the url's scheme.
It removes some of the boilerplate code required to download remote files from various providers:

.. code-block:: python

    from snakemake.remote import AUTO


    rule all:
        input:
            'foo'


    rule download:
        input:
            ftp_file_list=AUTO.remote([
                'ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxcat.tar.gz',
                'ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz'
            ], keep_local=True),
            http_file=AUTO.remote(
                'https://github.com/hetio/hetionet/raw/master/hetnet/tsv/hetionet-v1.0-nodes.tsv'
            )
        output:
            touch('foo')
        shell:
            """
            head {input.http_file}
            """
