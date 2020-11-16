# Executing this test case

To run this test, you need a running kubernetes setup.
For google cloud, see [here](https://snakemake.readthedocs.io/en/stable/executing/cloud.html#setup-kubernetes-on-google-cloud-engine).
With this, you can execute in case of google cloud:

    snakemake --kubernetes --use-conda --default-remote-provider GS --default-remote-prefix my-bucket

while replacing ``my-bucket`` with your storage bucket. The same test should also work on amazon (given that kubernetes is setup):

    snakemake --kubernetes --use-conda --default-remote-provider S3 --default-remote-prefix my-bucket
