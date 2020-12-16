.. _cloud:

===========================
Cloud execution
===========================



------------------------------------
Generic cloud support via Kubernetes
------------------------------------

Snakemake 4.0 and later supports execution in the cloud via Kubernetes.
This is independent of the cloud provider, but we provide the setup steps for GCE below.

Setup Kubernetes on Google cloud engine
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

First, install the `Google Cloud SDK <https://cloud.google.com/sdk/docs/quickstarts>`_.
Then, run

.. code-block:: console

    $ gcloud init

to setup your access.
Then, you can create a new kubernetes cluster via

.. code-block:: console

    $ gcloud container clusters create $CLUSTER_NAME --num-nodes=$NODES --scopes storage-rw

with ``$CLUSTER_NAME`` being the cluster name and ``$NODES`` being the number of cluster
nodes. If you intend to use google storage, make sure that ``--scopes storage-rw`` is set.
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
Importantly, this also means that you should not put large non-source files into the git repo, since Snakemake will try to upload them to kubernetes with every job.
With large files in the git repo, this can lead to performance issues or even random SSL errors from kubernetes.

It is further possible to forward arbitrary environment variables to the kubernetes
jobs via the flag ``--envvars`` (see ``snakemake --help``) or the ``envvars`` directive in the Snakefile.
The former should be used e.g. for platform specific variables (e.g. secrets that are only needed for your kubernetes setup), whereas the latter should be used for variables that are needed for the workflow itself, regardless of whether it is executed on kubernetes or with a different backend.

When executing, Snakemake will make use of the defined resources and threads
to schedule jobs to the correct nodes. In particular, it will forward memory requirements
defined as ``mem_mb`` to kubernetes. Further, it will propagate the number of threads
a job intends to use, such that kubernetes can allocate it to the correct cloud
computing node.


-------------------------------------------------------------
Executing a Snakemake workflow via Google Cloud Life Sciences
-------------------------------------------------------------

The `Google Cloud Life Sciences <https://cloud.google.com/life-sciences/docs/>`_
provides a rich application programming interface to design pipelines.
You'll first need to `follow instructions here <https://cloud.google.com/life-sciences/docs/quickstart>`_  to
create a Google Cloud Project and enable Life Sciences, Storage, and Compute Engine APIs,
and continue with the prompts to create credentials. You'll want to create
a service account for your host (it's easiest to give project Owner permissions), 
and save the json credentials. You'll want to export the full path to this file to ``GOOGLE_APPLICATION_CREDENTIALS`` :

.. code-block:: console

      $ export GOOGLE_APPLICATION_CREDENTIALS=$HOME/path/snakemake-credentials.json

If you lose the link to the credentials interface, you can `find it here <https://console.cloud.google.com/apis/credentials>`_.

Optionally, you can export ``GOOGLE_CLOUD_PROJECT`` as the name of your Google Cloud Project. By default, the project associated with your application credentials will be used.

.. code-block:: console

      $ export GOOGLE_CLOUD_PROJECT=my-project-name


Data in Google Storage
~~~~~~~~~~~~~~~~~~~~~~

Using this executor typically requires you to start with large data files
already in Google Storage, and then interact with them via the Google Storage
remote executor. An easy way to do this is to use the
`gsutil <https://cloud.google.com/storage/docs/uploading-objects>`_
command line client. For example, here is how we might upload a file
to storage using it:

.. code-block:: console

    $ gsutil -m cp mydata.txt gs://snakemake-bucket/1/mydata.txt

The ``-m`` parameter enables multipart uploads for large files, so you
can remove it if you are uploading one or more smaller files.
And note that you'll need to modify the file and bucket names.
Note that you can also easily use the Google Cloud Console interface, if
a graphical interface is preferable to you.

Environment Variables
~~~~~~~~~~~~~~~~~~~~~

**Important:** Google Cloud Life Sciences uses Google Compute, and does
**not** encrypt environment variables. If you specify environment
variables with the envvars directive or ``--envvars`` they will **not** be secrets.

Container Bases
~~~~~~~~~~~~~~~

By default, Google Life Sciences uses the latest stable version of
`snakemake/snakemake <https://hub.docker.com/r/snakemake/snakemake/tags>`_
on Docker Hub. You can choose to modify the container base with
the ``--container-image`` (or ``container_image`` from within Python),
however if you do so, your container must meet the following requirements:

 - have an entrypoint that can execute a ``/bin/bash`` command
 - have snakemake installed, either via ``conda activate snakemake`` or already on the path
 - also include snakemake Python dependencies for google.cloud

If you use any Snakemake container as a base, you should be good to go. If you'd
like to get a reference for requirements, it's helpful to look at the
`Dockerfile <https://github.com/snakemake/snakemake/blob/master/Dockerfile>`_
for Snakemake.

Requesting GPUs
~~~~~~~~~~~~~~~

The Google Life Sciences API currently has support for 
`NVIDIA GPUs <https://cloud.google.com/compute/docs/gpus#restrictions>`_, meaning that you can request a number of NVIDIA GPUs explicitly by adding ``nvidia_gpu`` or ``gpu`` to your Snakefile resources for a step:

.. code-block:: python

    rule a:
        output:
            "test.txt"
        resources:
            nvidia_gpu=1
        shell:
            "somecommand ..."

A specific `gpu model <https://cloud.google.com/compute/docs/gpus#introduction>`_ can be requested using ``gpu_model`` and lowercase identifiers like ``nvidia-tesla-p100`` or ``nvidia-tesla-p4``, for example: ``gpu_model="nvidia-tesla-p100"``. If you don't specify ``gpu`` or ``nvidia_gpu`` with a count, but you do specify a ``gpu_model``, the count will default to 1.

In addition to GPU for the Google Lifesciences Executor, you can request a `Google Cloud preemptible virtual machine <https://cloud.google.com/life-sciences/docs/reference/gcloud-examples#using_preemptible_vms>`_ for one or more steps. See the `rules documentation <https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#preemptible-virtual-machine>`_ for how to add one or more preemptible arguments.


Machine Types
~~~~~~~~~~~~~

To specify an exact `machine type <https://cloud.google.com/compute/docs/machine-types>`_
or a prefix to filter down to and then select based on other resource needs, 
you can set a default resource on the command line, either for a prefix or 
a full machine type:

.. code-block:: console

    --default-resources "machine_type=n1-standard"


If you want to specify the machine type as a resource, you can do that too:

.. code-block:: python

    rule a:
        output:
            "test.txt"
        resources:
            machine_type="n1-standard-8"
        shell:
            "somecommand ..."


If you request a gpu, this requires the "n1" prefix and your preference from
the file or command line will be overridden. Note that the default resources
for Google Life Sciences (memory and disk) are the same as for Tibanna.

Running the Life Sciences Executor
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When your Snakefile is ready, you can run snakemake to specify the life
sciences executor. Notice that we are also providing a remote prefix for our storage path,
along with a region.

.. code-block:: console

    $ snakemake --google-lifesciences --default-remote-prefix snakemake-testing-data --use-conda --google-lifesciences-region us-west1


For more details and examples, we recommend you reference the 
`Google Life Sciences Executor Tutorial <https://snakemake.readthedocs.io/en/stable/executor_tutorial/google_lifesciences.html>`_.


-----------------------------------------------------------------
Executing a Snakemake workflow via Tibanna on Amazon Web Services
-----------------------------------------------------------------

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

-----------------------------------------------------------------
Executing a Snakemake workflow via GA4GH TES
-----------------------------------------------------------------

The task execution service (`TES <https://github.com/ga4gh/task-execution-schemas>`_) is an application programming interface developed by the Global Alliance for Genomics and Health (`GA4GH <https://www.ga4gh.org/>`_).
It is used to process workflow tasks in a cloud environment.
A TES server can be easily implemented in a public cloud or at a commercial cloud provider.
Here, the TES standard provides an additional abstraction layer between the execution of a workflow (e.g. on your local machine) and technologies for execution of single tasks (e.g. based Kubernetes or HPC).
We recommend using either `Funnel <https://ohsu-comp-bio.github.io/funnel/>`_ or `TESK <https://github.com/EMBL-EBI-TSI/TESK/>`_  to install a TES server.
The guide here is based on Funnel (0.10.0).
To install and configure Funnel follow its official `documentation <https://ohsu-comp-bio.github.io/funnel/docs/>`_.

Configuration
~~~~~~~~~~~~~

Two steps are required to make a Snakemake workflow TES ready:

**Attach conda to rules:**
Execution of Snakemake tasks via TES means, Snakemake is running in a container in the cloud and it executes a specific rule (or a group of rules) with defined input/output data.
By default, the TES module uses the latest Snakemake container.
Running Snakemake within a container requires having all external tools installed within this container.
This can be done by providing a custom container image having installed Snakemake and other all required tools (e.g. BWA).
Or it can be done by attaching a conda environment to each rule, such that those tools will be installed within the running container.
For simplicity, this guide recommends to attach a specific conda environment to each rule, although it is more efficient in the long term to provide custom container images.

**Use remote files:**
The TES module requires using a remote file storage system for input/output files such that all files are available on the cloud machines and within their running container.
There are several options available in Snakemake to use remote files.
This guide recommends to use S3 (or SWIFT) object storage.

Execution
~~~~~~~~~

Funnel starts container in read only mode, which is good practice.
Anyhow, using the default Snakemake container image will likely require installing additional software within the running container.
Therefore, we need to set two conda specific variables such that new environments will be installed at `/tmp` which will be mounted as a writable volume in the container.

.. code-block:: console

    $ export CONDA_PKGS_DIRS=/tmp/conda
    $ export CONDA_ENVS_PATH=/tmp/conda

Next, using S3 or SWIFT storage, we also need to set credentials.

.. code-block:: console

    $ export AWS_ACCESS_KEY_ID=YOUR_ACCESS_KEY
    $ export AWS_SECRET_ACCESS_KEY=YOUR_SECRET_ACCESS_KEY

Now we can run Snakemake using:

.. code-block:: console

    $ snakemake \
        --tes $TES_URL \
        --use-conda \
        --envvars CONDA_PKGS_DIRS CONDA_ENVS_PATH AWS_ACCESS_KEY_ID AWS_SECRET_ACCESS_KEY \
        --conda-prefix $CONDA_ENVS_PATH \
        all
