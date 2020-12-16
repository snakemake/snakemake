.. _tutorial-azure-aks:

Auto-scaling Azure Kubernetes cluster without shared filesystem
---------------------------------------------------------------

.. _Snakemake: http://snakemake.readthedocs.io
.. _Python: https://www.python.org/

In this tutorial we will show how to execute a Snakemake workflow
on an auto-scaling Azure Kubernetes cluster without a shared file-system.
While Kubernetes is mainly known as microservice orchestration system with
self-healing properties, we will use it here simply as auto-scaling
compute orchestrator. One could use `persistent volumes in
Kubernetes <https://docs.microsoft.com/en-us/azure/aks/azure-files-dynamic-pv>`__
as shared file system, but this adds an unnecessary level of complexity
and most importantly costs. Instead we use cheap Azure Blob storage,
which is used by Snakemake to automatically stage data in and out for
every job.

Following the steps below you will

#. set up Azure Blob storage, download the Snakemake tutorial data and upload to Azure
#. then create an Azure Kubernetes (AKS) cluster
#. and finally run the analysis with Snakemake on the cluster 


Setup
:::::

To go through this tutorial, you need the following software installed:

* Python_ ≥3.5
* Snakemake_ ≥5.17

You should install conda as outlined in the :ref:`tutorial <tutorial-setup>`,
and then install full snakemake with:

.. code:: console

    conda create -c bioconda -c conda-forge -n snakemake snakemake

Make sure that the ``kubernetes`` and ``azure-storage-blob`` modules are installed
in this environment. Should they be missing install with:

.. code:: console

   pip install kubernetes
   pip install azure-storage-blob

In addition you will need the
`Azure CLI command <https://docs.microsoft.com/en-us/cli/azure/install-azure-cli?view=azure-cli-latest>`__ 
installed.

Create an Azure storage account and upload data
:::::::::::::::::::::::::::::::::::::::::::::::

We will be starting from scratch, i.e. we will 
create a new resource group and storage account. You can obviously reuse 
existing resources instead.

.. code:: console

   # change the following names as required
   # azure region where to run:
   region=southeastasia
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
A  more efficient way of uploading would be to use
`azcopy <https://docs.microsoft.com/en-us/azure/storage/common/storage-use-azcopy-v10>`__.

Create an auto-scaling Kubernetes cluster
:::::::::::::::::::::::::::::::::::::::::

.. code:: console

   # change the cluster name as you like
   clustername=snakemaks-aks
   az aks create --resource-group $resgroup --name $clustername \
       --vm-set-type VirtualMachineScaleSets --load-balancer-sku standard --enable-cluster-autoscaler \
       --node-count 1 --min-count 1 --max-count 3 --node-vm-size Standard_D3_v2

There is a lot going on here, so let’s unpack it: this creates an
`auto-scaling Kubernetes
cluster <https://docs.microsoft.com/en-us/azure/aks/cluster-autoscaler>`__
(``--enable-cluster-autoscaler``) called ``$clustername`` (i.e. ``snakemaks-aks``), which starts
out with one node (``--node-count 1``) and has a maximum of three nodes
(``--min-count 1 --max-count 3``). For real world applications you will
want to increase the maximum count and also increase the VM size. You
could for example choose a large instance from the DSv2 series and add a
larger disk with (``--node-osdisk-size``) if needed. See `here for more
info on Linux VM
sizes <https://docs.microsoft.com/en-us/azure/virtual-machines/linux/sizes>`__.

Note, if you are creating the cluster in the Azure portal, click on the
ellipsis under node-pools to find the auto-scaling option.

Next, let’s fetch the credentials for this cluster, so that we can
actually interact with it.

.. code:: console

   az aks get-credentials --resource-group $resgroup --name $clustername
   # print basic cluster info
   kubectl cluster-info



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


Now, we will need to setup the credentials that allow the Kubernetes nodes to
read and write from blob storage. For the AzBlob storage provider in
Snakemake this is done through the environment variables
``AZ_BLOB_ACCOUNT_URL`` and optionally ``AZ_BLOB_CREDENTIAL``. See the
`documentation <snakefiles/remote_files.html#microsoft-azure-storage>`__ for more info.
``AZ_BLOB_ACCOUNT_URL`` takes the form
``https://<accountname>.blob.core.windows.net`` and may also contain a
shared access signature (SAS), which is a powerful way to define fine grained
and even time controlled access to storage on Azure. The SAS can be part of the
URL, but if it’s missing, then you can set it with
``AZ_BLOB_CREDENTIAL`` or alternatively use the storage account key. To
keep things simple we’ll use the storage key here, since we already have it available,
but a SAS is generally more powerful. We’ll pass those variables on to the Kubernetes
with ``--envvars`` (see below).

Now you are ready to run the analysis:

.. code:: console

   export AZ_BLOB_ACCOUNT_URL="https://${stgacct}.blob.core.windows.net"
   export AZ_BLOB_CREDENTIAL="$stgkey"
   snakemake --kubernetes \
       --default-remote-prefix snakemake-tutorial --default-remote-provider AzBlob \
       --envvars AZ_BLOB_ACCOUNT_URL AZ_BLOB_CREDENTIAL --use-conda --jobs 3

This will use the default Snakemake image from Dockerhub. If you would like to use your
own, make sure that the image contains the same Snakemake version as installed locally
and also supports Azure Blob storage. If you plan to use your own image hosted on
 Azure Container Registries (ACR), make sure to attach the ACR to your Kubernetes 
 cluster. See `here <https://docs.microsoft.com/en-us/azure/aks/cluster-container-registry-integration>`__ for more info.


While Snakemake is running the workflow, it prints handy debug
statements per job, e.g.:

.. code:: console

   kubectl describe pod snakejob-c4d9bf9e-9076-576b-a1f9-736ec82afc64
   kubectl logs snakejob-c4d9bf9e-9076-576b-a1f9-736ec82afc64

With these you can also follow the scale-up of the cluster:

.. code:: console

   Events:
   Type     Reason             Age                From                Message
   ----     ------             ----               ----                -------
   Warning  FailedScheduling   60s (x3 over 62s)  default-scheduler   0/1 nodes are available: 1 Insufficient cpu.
   Normal   TriggeredScaleUp   50s                cluster-autoscaler  pod triggered scale-up: [{aks-nodepool1-17839284-vmss 1->3 (max: 3)}]

After a while you will see three nodes (each running one BWA job), which
was defined as the maximum above while creating your Kubernetes cluster:

.. code:: console

   $ kubectl get nodes
   NAME                                STATUS   ROLES   AGE   VERSION
   aks-nodepool1-17839284-vmss000000   Ready    agent   74m   v1.15.11
   aks-nodepool1-17839284-vmss000001   Ready    agent   11s   v1.15.11
   aks-nodepool1-17839284-vmss000002   Ready    agent   62s   v1.15.11

To get detailed information including historical data about used
resources, check Insights in the Azure portal under your AKS cluster
Monitoring/Insights. The alternative is an instant snapshot on the
command line:

::

   $ kubectl top node
   NAME                                CPU(cores)   CPU%   MEMORY(bytes)   MEMORY%
   aks-nodepool1-17839284-vmss000000   217m         5%     1796Mi          16%
   aks-nodepool1-17839284-vmss000001   1973m        51%    529Mi           4%
   aks-nodepool1-17839284-vmss000002   698m         18%    1485Mi          13%

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

Now that the execution is complete, the AKS cluster will scale down
automatically. If you are not planning to run anything else, it makes
sense to shut down it down entirely:

::

   az aks delete --name akscluster --resource-group $resgroup

