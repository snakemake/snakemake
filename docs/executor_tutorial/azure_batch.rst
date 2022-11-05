.. _tutorial-azure-batch:

Azure Batch Tutorial
---------------------------------------------------------------

.. _Snakemake: http://snakemake.readthedocs.io
.. _Python: https://www.python.org/

In this tutorial we will show how to execute a Snakemake workflow
on Azure batch nodes without a shared file-system. One could use attached storage 
solutions as a shared file system, but this adds an unnecessary level of complexity
and most importantly costs. Instead we use cheap Azure Blob storage,
which is used by Snakemake to automatically stage data in and out for
every job.

Following the steps below you will

#. set up Azure Blob storage, download the Snakemake tutorial data and upload to Azure
#. then create an Azure Batch account  
#. and finally run the analysis with Snakemake on the batch account


Setup
:::::

To go through this tutorial, you need the following software installed:

* Python_ ≥3.6
* Snakemake_ ≥7.18

You should install conda as outlined in the :ref:`tutorial <tutorial-setup>`,
and then install full snakemake with:

.. code:: console

    conda create -c bioconda -c conda-forge -n snakemake snakemake

Make sure that the ``azure-batch`` and ``azure-storage-blob`` modules are installed
in this environment. Should they be missing install with:

.. code:: console

   pip install azure-batch
   pip install azure-storage-blob

In addition you will need the
`Azure CLI command <https://docs.microsoft.com/en-us/cli/azure/install-azure-cli?view=azure-cli-latest>`__ 
installed.

Create an Azure storage account and upload example data
:::::::::::::::::::::::::::::::::::::::::::::::

We will be starting from scratch, i.e. we will 
create a new resource group and storage account. You can obviously reuse 
existing resources instead.

.. code:: console

   # change the following names as required
   # azure region where to run:
   region=westus

   # name of the resource group to create:
   resgroup=snakemaks-rg

   # name of storage account to create (all lowercase, no hyphens etc.):
   stgacct=snakemaksstg

   # create a resource group with name and in region as defined above
   az group create --name $resgroup --location $region

   # create a general purpose storage account with cheapest SKU
   az storage account create -n $stgacct -g $resgroup --sku Standard_LRS -l $region

Get a key for that account and save it as ``stgkey`` for later use:

.. code:: console

   stgkey=$(az storage account keys list -g $resgroup -n $storageacct | head -n1 | cut -f 3)

Next, you will create a storage container (think: bucket) to upload the Snakemake tutorial data to:

.. code:: console

   az storage container create --resource-group $resgroup --account-name $stgacct \
       --account-key $stgkey --name snakemake-tutorial

   cd /tmp

   git clone https://github.com/snakemake/snakemake-tutorial-data.git

   cd snakemake-tutorial-data

   az storage blob upload-batch -d snakemake-tutorial --account-name $stgacct \
       --account-key $stgkey -s data/ --destination-path data

We are using `az storage blob` for uploading, because that `az` is already installed.
Another cli tool for uploading to azure storage is 
`azcopy <https://docs.microsoft.com/en-us/azure/storage/common/storage-use-azcopy-v10>`__.

Create an Azure Batch Account
:::::::::::::::::::::::::::::::::::::::::

Create a new azure batch account. The batch account key will be given to snakemake to enable creation of batch resources by snakemake.

.. code:: console

    # can set variables as appropriate
    resgroup = snakemake-rg
    accountname = snakebatch01
    location = westus
    az batch account create --resource-group $resgroup --name $accountname --location $location\


.. code:: console

    az batch account keys --resource-group $resgroup --name $accountname --key-name primary



To run the test workflow, two primary environment variables need to be set local to the snakemake invocation. The azure batch account key, and the azure storage account url with an SAS key.

.. code:: console

     export AZ_BLOB_ACCOUNT_URL='${storage_account_url_with_sas}'
     export AZ_BATCH_ACCOUNT_KEY='${az_batch_account_key}'


Running the workflow
::::::::::::::::::::

Below we will task Snakemake to install software on the fly with conda.
For this we need a Snakefile with corresponding conda environment
yaml files. You can download the package containing all those files `here <workflow/snakedir.zip>`__.
After downloading, unzip it and cd into the newly created directory.

.. code:: console

   $ cd /tmp
   $ unzip ~/Downloads/snakedir.zip
   $ cd snakedir
   $ find .
   .
   ./Snakefile
   ./envs
   ./envs/calling.yaml
   ./envs/mapping.yaml


Now, we will need to setup the credentials that allow the Batch nodes to
read and write from blob storage. For the AzBlob storage provider in
Snakemake this is done through the environment variables
``AZ_BLOB_ACCOUNT_URL`` and optionally ``AZ_BLOB_CREDENTIAL``. See the
`documentation <snakefiles/remote_files.html#microsoft-azure-storage>`__ for more info.
``AZ_BLOB_ACCOUNT_URL`` takes the form ``https://<accountname>.blob.core.windows.net/`` 
or may also contain a shared access signature (SAS) ``https://<accountname>.blob.core.windows.net/<sas>``, 
which is a powerful way to define fine grained and even time controlled access to storage on Azure. 
The SAS can be part of the URL, but if it’s missing, then you can set it with
``AZ_BLOB_CREDENTIAL`` or alternatively use the storage account key. 
The SAS is generally a more powerful, and simple solution. We’ll pass those variables on to the batch nodes  
with ``--envvars`` (see below).

Now you are ready to run the analysis:

.. code:: console

    export AZ_BLOB_ACCOUNT_URL='${account_url_with_sas}'
    export AZ_BATCH_ACCOUNT_KEY='${az_batch_account_key}'

    snakemake \
     --jobs 3 \
     -rpf --default-remote-prefix snakemake-tutorial \
     --use-conda \
     --default-remote-provider AzBlob \
     --envvars AZ_BLOB_ACCOUNT_URL \
     --az-batch \
     --container-image snakemake/snakemake \
     --az-batch-account-url ${batch_account_url}


This will use the default Snakemake image from Dockerhub. If you would like to use your
own, make sure that the image contains the same Snakemake version as installed locally
and also supports Azure Blob storage. 

After completion all results including
logs can be found in the blob container. You will also find results
listed in the first Snakefile target downloaded to the working directoy.

::

   $ find snakemake-tutorial/
   snakemake-tutorial/
   snakemake-tutorial/calls
   snakemake-tutorial/calls/all.vcf


   $ az storage blob list  --container-name snakemake-tutorial --account-name $stgacct --account-key $stgkey -o table
   Name                     Blob Type    Blob Tier    Length    Content Type                       Last Modified              Snapshot
   -----------------------  -----------  -----------  --------  ---------------------------------  -------------------------  ----------
   calls/all.vcf            BlockBlob    Hot          90986     application/octet-stream           2020-06-08T05:11:31+00:00
   data/genome.fa           BlockBlob    Hot          234112    application/octet-stream           2020-06-08T03:26:54+00:00
   # etc.
   logs/mapped_reads/A.log  BlockBlob    Hot          346       application/octet-stream           2020-06-08T04:59:50+00:00
   mapped_reads/A.bam       BlockBlob    Hot          2258058   application/octet-stream           2020-06-08T04:59:50+00:00
   sorted_reads/A.bam       BlockBlob    Hot          2244660   application/octet-stream           2020-06-08T05:03:41+00:00
   sorted_reads/A.bam.bai   BlockBlob    Hot          344       application/octet-stream           2020-06-08T05:06:25+00:00
   # same for samples B and C

Once the execution is complete, the Batch nodes will scale down
automatically. If you are not planning to run anything else, it makes
sense to shut down it down entirely:

::

   az batch account delete --name $accountname --resource-group $resgroup