.. _tutorial-azure-batch:

Azure Batch Tutorial
---------------------------------------------------------------

.. _Snakemake: http://snakemake.readthedocs.io
.. _Python: https://www.python.org/
.. _AZCLI: https://docs.microsoft.com/en-us/cli/azure/install-azure-cli?view=azure-cli-latest

In this tutorial we will show how to execute a Snakemake workflow
on Azure Batch nodes using Azure Blob Storage. One could use attached storage 
solutions as a shared file system, but this adds an unnecessary level of complexity
and most importantly costs. Instead we use cheap Azure Blob Storage,
which is used by Snakemake to automatically stage data in and out for
every job. Please visit the `Azure Batch Documentation 
<https://learn.microsoft.com/en-us/azure/batch/batch-technical-overview#how-it-works>`__
for an overview of the various components of Azure Batch.

Following the steps below you will:

#. Set up Azure Blob Storage, and sync the Snakemake tutorial data to the storage container
#. Create an Azure Batch account  
#. Configure credentials
#. Run the example Sankemake workflow on the batch account


Setup
:::::

To go through this tutorial, you need the following software installed:

* Python_ ≥3.6
* Snakemake_ ≥7.28.0
* AZCLI_


First install conda as outlined in the :ref:`tutorial <tutorial-setup>`,
and then install the full Snakemake with some additional Azure related dependencies and AZCLI_:

.. code:: console

    conda create -c bioconda -c conda-forge -n snakemake snakemake msrest azure-batch azure-storage-blob azure-mgmt-batch azure-identity

Naturally, you can omit the deployment of such an environment in case you already have it, or you can update an existing Snakemake environment with the additional dependencies.

Create an Azure Storage Account and upload example data
:::::::::::::::::::::::::::::::::::::::::::::::

We will be starting from scratch, i.e. we will 
create a new resource group and storage account. You can obviously reuse 
existing resources instead.

.. code:: console

   # Login into Azure cli
   az login

   # change the following names as required
   # azure region where to run:
   export region=westus

   # name of the resource group to create:
   export resgroup=snakemaks-rg

   # name of storage account to create (all lowercase, no hyphens etc.):
   export stgacct=snakemaksstg

   # create a resource group with name and in region as defined above
   az group create --name $resgroup --location $region

   # create a general purpose storage account with cheapest SKU
   az storage account create -n $stgacct -g $resgroup --sku Standard_LRS -l $region

Get a key for that account and save it as ``stgkey``, then generate the storage account SAS token that expires one day later, you will use the SAS to authenticate to Blob storage:

.. code:: console

   # Get date 5 days from today
   export expiry_date=`date -u -d "+5 days" '+%Y-%m-%dT%H:%MZ'`

   # get the storage account key and storage endpoint
   export stgkey=$(az storage account keys list -g $resgroup -n $stgacct -o tsv | head -n1 | cut -f 4)
   export stgurl=$(az storage account show-connection-string -g $resgroup -n $stgacct --protocol https -o tsv | cut -f5,9 -d ';' | cut -f 2 -d '=')

   # get a storage account SAS token to use for AZ_BLOB_ACCOUNT_URL
   export sas=$(az storage account generate-sas --account-name $stgacct \
      --account-key $stgkey \
      --expiry $expiry_date \
      --https-only \
      --permissions acdlrw \
      --resource-types sco \
      --services bf
      --out tsv)

   # construct a blob account url with SAS token
   export storage_account_url_with_sas="${stgurl}?${sas}"

Next, you will create a storage container (think: bucket) to upload the Snakemake tutorial data to:

.. code:: console

   az storage container create --resource-group $resgroup --account-name $stgacct \
       --account-key $stgkey --name snakemake-tutorial

   cd /tmp

   git clone https://github.com/snakemake/snakemake-tutorial-data.git

   cd snakemake-tutorial-data

   az storage blob upload-batch -d snakemake-tutorial --account-name $stgacct \
       --account-key $stgkey -s data/ --destination-path data

Here we are using `az storage blob` for uploading the tutorial data, because the AZCLI_ is already installed.
Another cli tool for uploading to azure storage is 
`azcopy <https://docs.microsoft.com/en-us/azure/storage/common/storage-use-azcopy-v10>`__.

Azure Blob Storage Warning: 
:::::::::::::::::::

The snakemake azbatch executor will not work with data in a storage account that has "hierarchical namespace" enabled. 
Azure hierarchical namespace is a new api on azure storage that is also called "ADLS Gen2". 
Snakemake does not currently support this storage format because the Python API is distinct from traditional blob storage.
For more details see: https://learn.microsoft.com/en-us/azure/storage/blobs/data-lake-storage-namespace.


Create an Azure Batch Account
:::::::::::::::::::::::::::::::::::::::::

Create a new azure batch account and capture the batch account url and batch account key as environment variables. The batch account key will be given to snakemake to enable creation of batch resources by snakemake.

.. code:: console

    # can set variables as appropriate
    export accountname=snakebatch01
    az batch account create --resource-group $resgroup --name $accountname --location $region


The format of the batch account url is :code:`https://${accountname}.${region}.batch.azure.com`, which can be constructed from the output of the command :code:`az batch account list` or copied from the azure portal overview page of your batch account.

.. code:: console

    # get batch account url from command line
    export batch_endpoint=$(az batch account show --name $accountname --resource-group $resgroup --query "accountEndpoint" --output tsv)
    export batch_account_url="https://${batch_endpoint}"

    # set the batch account key
    export az_batch_account_key=$(az batch account keys list --resource-group $resgroup --name $accountname -o tsv | head -n1 | cut -f2)



To run the test workflow, two primary environment variables need to be set local to the snakemake invocation.
The azure batch account key, and the azure storage account url with an SAS credential. More details about the AZ_BLOB_ACCOUNT_URL 
are described in the section below. 

.. code:: console

     export AZ_BLOB_ACCOUNT_URL="${storage_account_url_with_sas}"
     export AZ_BATCH_ACCOUNT_KEY="${az_batch_account_key}"


Running the workflow
::::::::::::::::::::

Below we will run an example Snakemake workflow, using conda envrionments to install dependencies at runtime.
Clone the example workflow and cd into the directory:

.. code:: console

   $ git clone https://github.com/jakevc/snakemake-azbatch-example.git
   $ cd snakemake-azbatch-example
   $ tree 
   tree
    .
    ├── README.md
    ├── Snakefile
    ├── envs
    │   ├── calling.yaml
    │   ├── mapping.yaml
    │   └── stats.yaml
    ├── run.sh
    └── src
        └── plot-quals.py

To authenticate Azure Blob Storage, we set ``AZ_BLOB_ACCOUNT_URL`` 
which takes the form: ``https://<accountname>.blob.core.windows.net/?<sas_token>``. 
The SAS url can be constructed manually from the Azure portal, or on the command line using the commands shown in the above 
section on storage account configuration. The value for ``AZ_BLOB_ACCOUNT_URL`` must be enclosed in double quotes, as the SAS token 
contains special characters that need to be escaped.

When using azure storage and snakemake without the Azure Batch executor, it is valid to use storage account key credentials and the variable ``AZ_BLOB_CREDENTIAL``, 
but this type of authentication is not supported with Azure Batch so we must use ``AZ_BLOB_ACCOUNT_URL`` with an SAS token credential when using the Azure Batch executor.

We’ll pass the ``AZ_BLOB_ACCOUNT_URL`` on to the batch nodes with ``--envvars`` flag (see below). 

The following optional environment variables can be set to override their associated default values, 
and are used to change the runtime configuration of the batch nodes themselves:


.. list-table:: Optional Batch Node Configuration Environment Variables
   :widths: 40 40 40
   :header-rows: 1

   * - Environment Variable
     - Default Value
     - Description
   * - BATCH_POOL_IMAGE_PUBLISHER
     - microsoft-azure-batch
     - publisher of the vm image for the batch nodes 
   * - BATCH_POOL_IMAGE_OFFER
     - ubuntu-server-container
     - vm image offer for the batch nodes
   * - BATCH_POOL_IMAGE_SKU
     - 20-04-lts
     - vm image sku for batch nodes
   * - BATCH_POOL_VM_CONTAINER_IMAGE
     - ubuntu
     - batch nodes vm container image
   * - BATCH_POOL_VM_NODE_AGENT_SKU_ID
     - batch.node.ubuntu 20.04
     - sku id for batch node vm images
   * - BATCH_POOL_VM_SIZE
     - Standard_D2_v3
     - batch node vm image size
   * - BATCH_POOL_SUBNET_ID
     - None
     - subnetwork to deploy batch nodes into, requires the configuration of BATCH_MANAGED_IDENTITY
   * - BATCH_POOL_NODE_COUNT
     - 1
     - batch pool node count
   * - BATCH_POOL_RESOURCE_FILE_PREFIX
     - resource-files
     - container prefix for temporary resource files tar ball (Snakefile, envs)
   * - BATCH_NODE_START_TASK_SAS_URL
     - None
     - specify an SAS url to a bash script start task to run on each batch node
   * - BATCH_NODE_FILL_TYPE
     - spread
     - possible values ("spread", or "pack") 
   * - BATCH_NODE_COMMUNICATION_SIMPLIFIED 
     - None, "classic" 
     - If set, configures the batch pool to use the 'simplified' node communication mode. 
   * - BATCH_TASKS_PER_NODE
     - 1
     - the number of tasks allowed per batch node
   * - BATCH_MANAGED_IDENTITY_RESOURCE_ID
     - None
     - The resource ID of the managed identity to use
   * - BATCH_MANAGED_IDENTITY_CLIENT_ID
     - None
     - The client ID of the managed identity to use
   * - BATCH_CONTAINER_REGISTRY_URL
     - None
     - Container registry url to configure on the batch nodes 
   * - BATCH_CONTAINER_REGISTRY_USER
     - None
     - Container registry user, overrides managed identity authentication if set with password.
   * - BATCH_CONTAINER_REGISTRY_PASS
     - None
     - Container registry password
  
   

Now you are ready to run the analysis:

.. code:: console

    # required env variables
    export AZ_BLOB_PREFIX=snakemake-tutorial
    export AZ_BATCH_ACCOUNT_URL="${batch_account_url}"
    export AZ_BATCH_ACCOUNT_KEY="${az_batch_account_key}"
    export AZ_BLOB_ACCOUNT_URL="${storage_account_url_with_sas}"

    # optional environment variables with defaults listed

    # network and identity
    # export BATCH_POOL_SUBNET_ID=
    # export BATCH_MANAGED_IDENTITY_RESOURCE_ID=
    # export BATCH_MANAGED_IDENTITY_CLIENT_ID=

    # if unset, default is "classic"
    # export BATCH_NODE_COMMUNICATION_SIMPLIFIED=true

    # don't recommend changing 
    # export BATCH_POOL_IMAGE_PUBLISHER=microsoft-azure-batch
    # export BATCH_POOL_IMAGE_OFFER=ubuntu-server-container
    # export BATCH_POOL_IMAGE_SKU=20-04-lts
    # export BATCH_POOL_RESOURCE_FILE_PREFIX=resource-files

    # export BATCH_POOL_VM_CONTAINER_IMAGE=ubuntu
    # export BATCH_POOL_VM_NODE_AGENT_SKU_ID="batch.node.ubuntu 20.04"

    # can be used to add a startup task to the batch nodes formatted as an sas url to a bash script
    # export BATCH_NODE_START_TASK_SAS_URL=

    # can be useful to alter task distribution across nodes

    # export BATCH_POOL_VM_SIZE=Standard_D2_v3
    # export BATCH_NODE_FILL_TYPE=spread
    # export BATCH_POOL_NODE_COUNT=1
    # export BATCH_TASKS_PER_NODE=1

    # container registry configuration to pull container image from custom registry
    # export BATCH_CONTAINER_REGISTRY_URL=
    # export BATCH_CONTAINER_REGISTRY_USER=
    # export BATCH_CONTAINER_REGISTRY_PASS=

    snakemake \
        --jobs 3 \
        -rpf --verbose --default-remote-prefix $AZ_BLOB_PREFIX \
        --use-conda \
        --default-remote-provider AzBlob \
        --envvars AZ_BLOB_ACCOUNT_URL \
        --az-batch \
        --container-image snakemake/snakemake \
        --az-batch-account-url $AZ_BATCH_ACCOUNT_URL

This will use the default Snakemake image from Dockerhub. If you would like to use your
own, make sure that the image contains the same Snakemake version as installed locally
and also supports Azure Blob Storage. The optional BATCH_CONTAINER_REGISTRY can be configured 
to fetch from your own container registry. If that registry is an Azure Container Registry 
that the managed identity has access to, then the BATCH_CONTAINER_REGISTRY_USER and BATCH_CONTAINER_REGISTRY_PASS is not needed. 

After completion all results including
logs can be found in the blob container prefix specified by `--default-remote-prefix`.

::

   $ az storage blob list  --container-name snakemake-tutorial --account-name $stgacct --account-key $stgkey -o table
   Name                                                                                            IsDirectory    Blob Type    Blob Tier    Length    Content Type              Last Modified              Snapshot
  ----------------------------------------------------------------------------------------------  -------------  -----------  -----------  --------  ------------------------  -------------------------  ----------
  data/genome.fa                                                                                                 BlockBlob    Hot          234112    application/octet-stream  2022-12-14T23:28:00+00:00
  data/genome.fa.amb                                                                                             BlockBlob    Hot          2598      application/octet-stream  2022-12-14T23:28:01+00:00
  data/genome.fa.ann                                                                                             BlockBlob    Hot          83        application/octet-stream  2022-12-14T23:28:01+00:00
  data/genome.fa.bwt                                                                                             BlockBlob    Hot          230320    application/octet-stream  2022-12-14T23:28:01+00:00
  data/genome.fa.fai                                                                                             BlockBlob    Hot          18        application/octet-stream  2022-12-14T23:28:01+00:00
  data/genome.fa.pac                                                                                             BlockBlob    Hot          57556     application/octet-stream  2022-12-14T23:28:00+00:00
  data/genome.fa.sa                                                                                              BlockBlob    Hot          115160    application/octet-stream  2022-12-14T23:28:01+00:00
  data/samples/A.fastq                                                                                           BlockBlob    Hot          5752788   application/octet-stream  2022-12-14T23:28:04+00:00
  data/samples/B.fastq                                                                                           BlockBlob    Hot          5775000   application/octet-stream  2022-12-14T23:28:06+00:00
  data/samples/C.fastq                                                                                           BlockBlob    Hot          5775000   application/octet-stream  2022-12-14T23:28:02+00:00
  logs/mapped_reads/A.log                                                                                        BlockBlob    Hot                    application/octet-stream  2022-12-28T18:14:33+00:00
  logs/mapped_reads/B.log                                                                                        BlockBlob    Hot                    application/octet-stream  2022-12-28T18:15:25+00:00
  logs/mapped_reads/C.log                                                                                        BlockBlob    Hot                    application/octet-stream  2022-12-28T18:16:17+00:00
  results/calls/all.vcf                                                                                          BlockBlob    Hot          90962     application/octet-stream  2022-12-28T18:22:20+00:00
  results/mapped_reads/A.bam                                                                                     BlockBlob    Hot          2258050   application/octet-stream  2022-12-28T18:14:33+00:00
  results/mapped_reads/B.bam                                                                                     BlockBlob    Hot          2262766   application/octet-stream  2022-12-28T18:15:25+00:00
  results/mapped_reads/C.bam                                                                                     BlockBlob    Hot          2262766   application/octet-stream  2022-12-28T18:16:17+00:00
  results/plots/quals.svg                                                                                        BlockBlob    Hot          12571     application/octet-stream  2022-12-28T19:16:28+00:00
  results/sorted_reads/A.bam                                                                                     BlockBlob    Hot          2244652   application/octet-stream  2022-12-28T18:17:10+00:00
  results/sorted_reads/A.bam.bai                                                                                 BlockBlob    Hot          344       application/octet-stream  2022-12-28T18:19:48+00:00
  results/sorted_reads/B.bam                                                                                     BlockBlob    Hot          2248758   application/octet-stream  2022-12-28T18:18:08+00:00
  results/sorted_reads/B.bam.bai                                                                                 BlockBlob    Hot          344       application/octet-stream  2022-12-28T18:20:36+00:00
  results/sorted_reads/C.bam                                                                                     BlockBlob    Hot          2248758   application/octet-stream  2022-12-28T18:18:58+00:00
  results/sorted_reads/C.bam.bai                                                                                 BlockBlob    Hot          344       application/octet-stream  2022-12-28T18:21:23+00:00

Once the execution is complete, the batch nodes will scale down
automatically. If you are not planning to run anything else, it makes
sense to shut down it down entirely:

::

   az batch account delete --name $accountname --resource-group $resgroup


Defining a Start Task
:::::
A start task can be optionally specified as a shell scirpt that runs during each node's startup as it's added to the batch pool.
To specify a start task, set the environment variable BATCH_NODE_START_TASK_SAS_URL to the SAS url of a start task shell script.
Store your shell script in a blob storage account and generate an SAS url to a shell script blob object. 
You can generate an SAS URL to the blob using the azure portal or the command line using the following command structure: 

::

  container="container-name"
  expiry="2024-01-01"
  blob_name="starttask.sh"
  SAS_TOKEN=$(az storage blob generate-sas --account-name $stgacct --container-name $container --name $blob_name --permissions r --auth-mode login --as-user --expiry $expiry -o tsv)
  BLOB_URL=$(az storage blob url --account-name cromwellstorage --container-name snaketest --name starttask.sh --auth-mode login -o tsv)

  # then export the full SAS URL
  export BATCH_NODE_START_TASK_SAS_URL="${BLOB_URL}?${SAS_TOKEN}"


Autoscaling and Task Distribution
:::::

The azure batch executor supports autoscaling of the batch nodes by including the flag ``--az-batch-enable-autoscale``. 
This flag sets the initial dedicated node count of the pool to zero, and re-evaluates the number of nodes to be spun up or down based on the number of remaining tasks to run over a five minute interval. 
Since five minutes is the smallest allowed interval for azure batch autoscaling, this feature becomes more useful for long running jobs. For more information on azure batch autoscaling configuration, see: https://learn.microsoft.com/en-us/azure/batch/batch-automatic-scaling.

For shorter running jobs it might be more cost/time effective to set VM size with more cores (`BATCH_POOL_VM_SIZE`) and increase the number of `BATCH_TASKS_PER_NODE`. Or, if you want to keep tasks running on separate nodes, you can set a larger number for `BATCH_POOL_NODE_COUNT`. 
It may require experimentation to find the most efficient/cost effective task distribution model for your use case depending on what you are optimizing for. For more details on limitations of azure batch node / task distribution see: https://learn.microsoft.com/en-us/azure/batch/batch-parallel-node-tasks.