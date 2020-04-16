
.. _tutorial-google-lifesciences:

Google Life Sciences Tutorial
------------------------------

.. _Snakemake: http://snakemake.readthedocs.io
.. _Snakemake Remotes: https://snakemake.readthedocs.io/en/stable/snakefiles/remote_files.html


Setup
:::::

To go through this tutorial, you need the following software installed:

* Python_ ≥3.5
* Snakemake_ 5.2.3

You can follow the instructions on the :ref:`setup <tutorial-setup>` page,
and return here when you are in a new virtual environment with snakemake installed.

Install Dependencies
::::::::::::::::::::

You can now install dependencies for the Google Life Sciences Executor. This
is represented fully with:

.. code:: console

    pip install snakemake[google-cloud]


which should include all of the following:

.. code:: console

    pip install --upgrade google-api-python-client
    pip install --upgrade google-cloud-storage
    pip install oauth2client
    pip install crc32c


Credentials
:::::::::::

Snakemake requires `GOOGLE_APPLICATION_CREDENTIALS`, and since you might want to
run this is (non Google places) too, you should `download your service account <https://console.cloud.google.com/iam-admin/iam>`_
key and export it to the environment.

.. code:: console

    export GOOGLE_APPLICATION_CREDENTIALS="/home/[username]/credentials.json"



Step 1: Upload Your Data
::::::::::::::::::::::::

We will be obtaining inputs from Google Cloud Storage, as well as saving
outputs there. You should first clone the repository with the Snakemake tutorial data:


.. code:: console

    git clone https://github.com/snakemake/snakemake-tutorial-data
    cd snakemake-tutorial-data


And then either manually create a bucket and upload data files there, or
use the `provided script and instructions <https://github.com/snakemake/snakemake-tutorial-data#google-cloud-storage>`_
do do it programatically from the command line. As an example:


.. code:: console

    export GOOGLE_APPLICATION_CREDENTIALS="/path/to/credentials.json"
    python upload_google_storage.py my-snakemake-bucket/data  


Your bucket (and the folder prefix) will be referred to as the
`--default-remote-prefix` when you run snakemake.


Step 2: Write your Snakefile, Environment File, and Scripts
:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

Now that we've exported our credentials and have all dependencies installed, let's
get our workflow! This is the exact same workflow from the :ref:`basic tutorial<tutorial-basics>`,
so if you need a refresher on the design or basics, please see those pages.
We won't need to clone data from `snakemake-tutorial-data <https://github.com/snakemake/snakemake-tutorial-data>`_
because we will be using publicly accessible data on Google Storage.

First, how does a working directory work for this executor? The present
working directory, as identified by Snakemake that has the Snakefile, and where
a more advanced setup might have a folder of environment specifications (env) a folder of scripts 
(scripts), and rules (rules), is considered within the context of the build.
When the Google Life Sciences executor is used, it generates a build package of all
of the files here (within a reasonable size) and uploads those to storage. This
package includes the .snakemake folder that would have been generated locally.
The build package is then downloaded and extracted by each cloud executor, which
is a Google Compute instance.

We next need an `environment.yml` file that will define the dependencies
that we want installed with conda for our job. If you cloned the "snakemake-tutorial-data"
repository you will already have this, and you are good to go. If not, save this to `environment.yml`
in your working directory:

.. code:: yaml

    channels:
      - conda-forge
      - bioconda
    dependencies:
      - bioconda::snakemake-minimal =5.4.5
      - python =3.6
      - jinja2 =2.10
      - networkx =2.1
      - matplotlib =2.2.3
      - graphviz =2.38.0
      - bcftools =1.9
      - samtools =1.9
      - bwa =0.7.17
      - pysam =0.15.0
    

Notice that we reference this `environment.yml` file in the Snakefile below.
Importantly, if you were optimizing a pipeline, you would likely have a folder
"envs" with more than one environment specification, one for each step.
This workflow uses the same environment (with many dependencies) instead of
this strategy to minimize the number of files for you.

The Snakefile then has the following content. It's important to note
that we have not customized this file from the basic tutorial to hard code 
any storage or executor. We will be telling snakemake to change the executor 
via command line arguments.

.. code:: python

    SAMPLES = ["A", "B"]

    rule all:
        input:
            "plots/quals.svg"

    rule bwa_map:
        input:
            "genome.fa",
            "samples/{sample}.fastq"
        conda:
            "environment.yml"
        output:
            "mapped_reads/{sample}.bam"
        shell:
            "bwa mem {input} | samtools view -Sb - > {output}"

    rule samtools_sort:
        input:
            "mapped_reads/{sample}.bam"
        output:
            "sorted_reads/{sample}.bam"
        conda:
            "environment.yml"
        shell:
            "samtools sort -T sorted_reads/{wildcards.sample} "
            "-O bam {input} > {output}"

    rule samtools_index:
        input:
            "sorted_reads/{sample}.bam"
        output:
            "sorted_reads/{sample}.bam.bai"
        conda:
            "environment.yml"
        shell:
            "samtools index {input}"

    rule bcftools_call:
        input:
            fa="genome.fa",
            bam=expand("sorted_reads/{sample}.bam", sample=SAMPLES),
            bai=expand("sorted_reads/{sample}.bam.bai", sample=SAMPLES)
        output:
            "calls/all.vcf"
        conda:
            "environment.yml"
        shell:
            "samtools mpileup -g -f {input.fa} {input.bam} | "
            "bcftools call -mv - > {output}"

    rule plot_quals:
        input:
            "calls/all.vcf"
        output:
            "plots/quals.svg"
        conda:
            "environment.yml"
        script:
            "plot-quals.py"



And let's also write the script in our present working directory for the last step
to do the plotting - call this `plot-quals.py`:

.. code:: python

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from pysam import VariantFile

    quals = [record.qual for record in VariantFile(snakemake.input[0])]
    plt.hist(quals)

    plt.savefig(snakemake.output[0])


Step 3: Run Snakemake
:::::::::::::::::::::

Now let's run Snakemake with the Google Life Sciences Executor.


.. code:: console

    snakemake --google-lifesciences --default-remote-prefix snakemake-testing-data --use-conda --google-lifesciences-region us-west1 --container-image snakemake/snakemake:v5.10.0


The flags above refer to:

 - `--google-lifesciences`: to indicate that we want to use the Google Life Sciences API
 - `--default-remote-prefix`: refers to the Google Storage bucket. The bucket name is "snakemake-testing-data" and the "subfolder" (or path) (not defined above) would be a subfolder, if needed.
 - `--google-lifesciences-region`: the region that you want the instances to deploy to. Your storage bucket should be accessible from here, and your selection can have a small influence on the machine type selected.


Once you submit the job, you'll immediately see the familiar Snakemake console output,
but with additional lines for inspecting google compute instances with gcloud:

.. code:: console

    Building DAG of jobs...
    Unable to retrieve additional files from git. This is not a git repository.
    Using shell: /bin/bash
    Rules claiming more threads will be scaled down.
    Job counts:
    	count	jobs
    	1	all
    	1	bcftools_call
    	2	bwa_map
	1	plot_quals
	2	samtools_index
	2	samtools_sort
	9

    [Thu Apr 16 19:16:24 2020]
    rule bwa_map:
        input: snakemake-testing-data/genome.fa, snakemake-testing-data/samples/B.fastq
        output: snakemake-testing-data/mapped_reads/B.bam
        jobid: 8
        wildcards: sample=B
        resources: mem_mb=15360, disk_mb=128000

    Get status with:
    gcloud config set project snakemake-testing
    gcloud beta lifesciences operations describe 13586583122112209762
    gcloud beta lifesciences operations list


Take not of those last three lines to describe and list operations - this is how you
get complete error and output logs for the run, which we will demonstrate using later.


Step 4: View Results
::::::::::::::::::::



Step 5: Debugging
:::::::::::::::::

Let's introduce an error (purposefully) into our Snakefile to practice debugging.
Let's remove the conda environment.yml file for the first rule, so we would
expect that Snakemake won't be able to find the executables for bwa and samtools.
In your Snakefile, change this:

.. code:: python

    rule bwa_map:
        input:
            "genome.fa",
            "samples/{sample}.fastq"
        conda:
            "environment.yml"
        output:
            "mapped_reads/{sample}.bam"
        shell:
            "bwa mem {input} | samtools view -Sb - > {output}"


to this:

.. code:: python

    rule bwa_map:
        input:
            "genome.fa",
            "samples/{sample}.fastq"
        output:
            "mapped_reads/{sample}.bam"
        shell:
            "bwa mem {input} | samtools view -Sb - > {output}"


And then run the job again:

.. code:: console

    snakemake --google-lifesciences --default-remote-prefix snakemake-testing-data --use-conda --google-lifesciences-region us-west1 --container-image snakemake/snakemake:v5.10.0
