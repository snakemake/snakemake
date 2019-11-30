===========================
Cluster and cloud execution
===========================

.. _cloud:

-------------
Cloud Support
-------------

Snakemake 4.0 and later supports execution in the cloud via Kubernetes.
This is independent of the cloud provider, but we provide the setup steps for GCE below.

Google cloud engine
~~~~~~~~~~~~~~~~~~~

First, install the `Google Cloud SDK <https://cloud.google.com/sdk/docs/quickstarts>`_.
Then, run

.. code-block:: console

    $ gcloud init

to setup your access.
Then, you can create a new kubernetes cluster via

.. code-block:: console

    $ gcloud container clusters create $CLUSTER_NAME --num-nodes=$NODES --scopes storage-rw

with ``$CLUSTER_NAME`` being the cluster name and ``$NODES`` being the number of cluster
nodes. If you intent to use google storage, make sure that `--scopes storage-rw` is set.
This enables Snakemake to write to the google storage from within the cloud nodes.
Next, you configure Kubernetes to use the new cluster via

.. code-block:: console

    $ gcloud container clusters get-credentials $CLUSTER_NAME


If you are having issues with authentication, please refer to the help text:

.. code-block:: console

    $ gcloud container clusters get-credentials --help

You likely also want to use google storage for reading and writing files.
For this, you will additionally need to authenticate with your google cloud account via

.. code-block:: console

    $ gcloud auth application-default login

This enables Snakemake to access google storage in order to check existence and modification dates of files.
Now, Snakemake is ready to use your cluster.

**Important:** After finishing your work, do not forget to delete the cluster with

.. code-block:: console

    $ gcloud container clusters delete $CLUSTER_NAME

in order to avoid unnecessary charges.


.. _kubernetes:


Executing a Snakemake workflow via kubernetes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Assuming that kubernetes has been properly configured (see above), you can
execute a workflow via:

.. code-block:: console

    snakemake --kubernetes --use-conda --default-remote-provider $REMOTE --default-remote-prefix $PREFIX

In this mode, Snakemake will assume all input and output files to be stored in a given
remote location, configured by setting ``$REMOTE`` to your provider of choice
(e.g. ``GS`` for Google cloud storage or ``S3`` for Amazon S3) and ``$PREFIX``
to a bucket name or subfolder within that remote storage.
After successful execution, you find your results in the specified remote storage.
Of course, if any input or output already defines a different remote location, the latter will be used instead.
Importantly, this means that Snakemake does **not** require a shared network
filesystem to work in the cloud.


.. sidebar:: Note

  Consider to :ref:`group jobs <snakefiles-grouping>` in order to minimize overhead, in particular for short-running jobs.

Currently, this mode requires that the Snakemake workflow is stored in a git repository.
Snakemake uses git to query necessary source files (the Snakefile, scripts, config, ...)
for workflow execution and encodes them into the kubernetes job.

It is further possible to forward arbitrary environment variables to the kubernetes
jobs via the flag ``--kubernetes-env`` (see ``snakemake --help``).

When executing, Snakemake will make use of the defined resources and threads
to schedule jobs to the correct nodes. In particular, it will forward memory requirements
defined as `mem_mb` to kubernetes. Further, it will propagate the number of threads
a job intends to use, such that kubernetes can allocate it to the correct cloud
computing node.


Executing a Snakemake workflow via Google Cloud Life Sciences
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The `Google Cloud Life Sciences <https://cloud.google.com/life-sciences/docs/>`_
provides a rich application programming interface to design pipelines.
You'll first need to `follow instructions here <https://cloud.google.com/life-sciences/docs/quickstart>`_  to
create a Google Cloud Project and enable Life Sciences, Storage, and Compute Engine APIs,
and continue with the prompts to create credentials. You'll want to create
a service account for your host (it's easiest to give project Owner permissions), 
and save the json credentials. You'll want to export the full path to this file to `GOOGLE_APPLICATION_CREDENTIALS`:

.. code-block:: console

      $ export GOOGLE_APPLICATION_CREDENTIALS=$HOME/path/snakemake-credentials.json

If you lose the link to the credentials interface, you can `find it here <https://console.cloud.google.com/apis/credentials>`_.

Optionally, you can export `GOOGLE_CLOUD_PROJECT` as the name of your Google
Cloud Project. By default, the project associated with your application 
credentials will be used.

.. code-block:: console

      $ export GOOGLE_CLOUD_PROJECT=my-project-name

The dependencies that you'll need for snakemake are:

 - gcc
 - python dev
 - google cloud python client libraries
 - oauth2client


Data in Google Storage
::::::::::::::::::::::

Using this executor typically requires you to start with large data files
already in Google Storage, and then interact with them via the Google Storage
remote executor. An easy way to do this is to use the
`gsutil <https://cloud.google.com/storage/docs/uploading-objects>`_
command line client. For example, here is how we might upload a file
to storage using it:


.. code-block:: console

      $ gsutil -m cp mydata.txt gs://snakemake-bucket/1/mydata.txt

The `-m` parameter enables multipart uploads for large files, so you
can remove it if you are uploading one or more smaller files.
And note that you'll need to modify the file and bucket names.
If you want a programmatic method in Python, here is an example script 
that shows you how you might want to upload files in advance using the 
Google API Python Client. 

.. code-block:: python

    from google.cloud import storage
    import glob
    import os

    client = storage.Client()
    here = os.path.dirname(os.path.abspath(__file__))

    bucket_name = "snakemake-testing"

    try:
        bucket = client.get_bucket(bucket_name)
    except:
        bucket = client.create_bucket(bucket_name)

    # key, full path, value: upload path in storage
    filenames = dict()

    for x in os.walk(os.path.join(here, 'scripts')):
        for name in glob.glob(os.path.join(x[0], '*')):
            if not os.path.isdir(name):
                filenames[name] = name.replace(here + os.path.sep, '')

    # Upload files in script and base gooogle life sciences folder
    for filename, remotename in filenames.items():
        blob = bucket.blob(remotename)

        if not blob.exists():
            print("Uploading %s to %s" %(remotename, bucket_name))
            blob.upload_from_filename(filename, content_type="text/plain")


Note that you can also easily use the Google Cloud Console interface, if
a graphical interface is preferable to you.

Environment Variables
:::::::::::::::::::::

**Important:** Google Cloud Life Sciences uses Google Compute, and does
**not** encrypt environment variables. If you specify environment
variables with `--envars` they will **not** be secrets.


Container Bases
:::::::::::::::

By default, Google Life Sciences uses the latest stable version of
`snakemake/snakemake <https://hub.docker.com/r/snakemake/snakemake/tags>`_
on Docker Hub. You can choose to modify the container base with
the `--container-image` (or `container_image` from within Python), 
however if you do so, your container must meet the following requirements:

 - have an entrypoint that can execute a `/bin/bash` command
 - have snakemake installed, either via `source conda snakemake` or already on the path
 - also include snakemake Python dependencies for google.cloud

If you use any Snakemake container as a base, you should be good to go. If you'd
like to get a reference for requirements, it's helpful to look at the
`Dockerfile <https://github.com/snakemake/snakemake/blob/master/Dockerfile>`_
for Snakemake.


Executing a Snakemake workflow via Tibanna on Amazon Web Services
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

First, install `Tibanna <https://tibanna.readthedocs.io/en/latest/>`_.

.. code-block:: console

    $ pip install -U tibanna


Set up aws configuration either by creating files ``~/.aws/credentials`` and ``~/.aws/config`` 
or by setting up environment variables as below (see Tibanna or AWS documentation for more details):

.. code-block:: console

    $ export AWS_ACCESS_KEY_ID=<AWS_ACCESS_KEY>
    $ export AWS_SECRET_ACCESS_KEY=<AWS_SECRET_ACCESS_KEY>
    $ export AWS_DEFAULT_REGION=<AWS_DEFAULT_REGION>


As an AWS admin, deploy Tibanna Unicorn to Cloud with permissions to a specific S3 bucket.
Name the Unicorn / Unicorn usergroup with the ``--usergroup`` option.
Unicorn is a serverless scheduler, and keeping unicorn on the cloud does not incur extra cost. 
One may have many different unicorns with different names and different bucket permissions.
Then, add other (IAM) users to the user group that has permission to use this unicorn / buckets.

.. code-block:: console

    $ tibanna deploy_unicorn -g <name> -b <bucket>
    $ tibanna add_user -u <username> -g <name>


As a user that has been added to the group (or as an admin), set up the default unicorn.

.. code-block:: console

    $ export TIBANNA_DEFAULT_STEP_FUNCTION_NAME=tibanna_unicorn_<name>


Then, you can run as many snakemake runs as you wish as below, inside a directory that contains
Snakefile and other necessary components (e.g. ``env.yml``, ``config.json``, ...).

.. code-block:: console

    $ snakemake --tibanna --default-remote-prefix=<bucketname>/<subdir> [<other options>]


In this mode, Snakemake will assume all input and output files to be stored in the specified remote location
(a subdirectory inside a given S3 bucket.)
After successful execution, you find your results in the specified remote storage.
Of course, if any input or output already defines a different remote location, the latter will be used instead.
In that case, Tibanna Unicorn must be deployed with all the relevant buckets (``-b bucket1,bucket2,bucket3,...``)
to allow access to the Unicorn serverless components.
Snakemake will assign 3x of the total input size as the allocated space for each execution. The execution may fail
if the total input + output + temp file sizes exceed this estimate.

In addition to regular snakemake options, ``--precommand=<command>`` option allows sending a command to execute before
executing each step on an isolated environment. This kind of command could involve downloading or installing
necessary files that cannot be handled using conda (e.g. the command may begin with ``wget``, ``git clone``, etc.) 


To check Tibanna execution logs, first use ``tibanna stat`` to see the list of all the individual runs.

.. code-block:: console

    $ tibanna stat -n <number_of_executions_to_view> -l


Then, check the detailed log for each job using the Tibanna job id that can be obtained from the first column
of the output of ``tibanna stat``.


.. code-block:: console

    $ tibanna log -j <jobid>


.. sidebar:: Note

  Consider to :ref:`group jobs <snakefiles-grouping>` in order to minimize overhead, in particular for short-running jobs.


When executing, Snakemake will make use of the defined resources and threads
to schedule jobs to the correct nodes. In particular, it will forward memory requirements
defined as `mem_mb` to Tibanna. Further, it will propagate the number of threads
a job intends to use, such that Tibanna can allocate it to the most cost-effective
cloud compute instance available.


.. _cluster:

-----------------
Cluster Execution
-----------------


Snakemake can make use of cluster engines that support shell scripts and have access to a common filesystem, (e.g. the Sun Grid Engine).
In this case, Snakemake simply needs to be given a submit command that accepts a shell script as first positional argument:

.. code-block:: console

    $ snakemake --cluster qsub -j 32


Here, ``-j`` denotes the number of jobs submitted being submitted to the cluster at the same time (here 32).
The cluster command can be decorated with job specific information, e.g.

.. sidebar:: Note

  Consider to :ref:`group jobs <snakefiles-grouping>` in order to minimize overhead, in particular for short-running jobs.


.. code-block:: console

    $ snakemake --cluster "qsub {threads}"

Thereby, all keywords of a rule are allowed (e.g. rulename, params, input, output, threads, priority, ...).
For example, you could encode the expected running time into params:

.. code-block:: python

    rule:
        input:  ...
        output: ...
        params: runtime="4h"
        shell: ...

and forward it to the cluster scheduler:

.. code-block:: console

    $ snakemake --cluster "qsub --runtime {params.runtime}"

If your cluster system supports `DRMAA <http://www.drmaa.org/>`_, Snakemake can make use of that to increase the control over jobs.
E.g. jobs can be cancelled upon pressing ``Ctrl+C``, which is not possible with the generic ``--cluster`` support.
With DRMAA, no ``qsub`` command needs to be provided, but system specific arguments can still be given as a string, e.g.

.. code-block:: console

    $ snakemake --drmaa " -q username" -j 32

Note that the string has to contain a leading whitespace.
Else, the arguments will be interpreted as part of the normal Snakemake arguments, and execution will fail.


Job Properties
~~~~~~~~~~~~~~

When executing a workflow on a cluster using the ``--cluster`` parameter (see below), Snakemake creates a job script for each job to execute. This script is then invoked using the provided cluster submission command (e.g. ``qsub``). Sometimes you want to provide a custom wrapper for the cluster submission command that decides about additional parameters. As this might be based on properties of the job, Snakemake stores the job properties (e.g. name, rulename, threads, input, output, params etc.) as JSON inside the job script (for group jobs, the rulename will be "GROUP", otherwise it will be the same as the job name). For convenience, there exists a parser function `snakemake.utils.read_job_properties` that can be used to access the properties. The following shows an example job submission wrapper:

.. code-block:: python

    #!python

    #!/usr/bin/env python3
    import os
    import sys

    from snakemake.utils import read_job_properties

    jobscript = sys.argv[1]
    job_properties = read_job_properties(jobscript)

    # do something useful with the threads
    threads = job_properties[threads]

    # access property defined in the cluster configuration file (Snakemake >=3.6.0)
    job_properties["cluster"]["time"]

    os.system("qsub -t {threads} {script}".format(threads=threads, script=jobscript))

.. _getting_started-visualization:

-------------
Visualization
-------------

To visualize the workflow, one can use the option ``--dag``.
This creates a representation of the DAG in the graphviz dot language which has to be postprocessed by the graphviz tool ``dot``.
E.g. to visualize the DAG that would be executed, you can issue:

.. code-block:: console

    $ snakemake --dag | dot | display

For saving this to a file, you can specify the desired format:

.. code-block:: console

    $ snakemake --dag | dot -Tpdf > dag.pdf

To visualize the whole DAG regardless of the eventual presence of files, the ``forceall`` option can be used:

.. code-block:: console

    $ snakemake --forceall --dag | dot -Tpdf > dag.pdf

Of course the visual appearance can be modified by providing further command line arguments to ``dot``.
