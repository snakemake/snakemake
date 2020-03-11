# Instruction for testing of Azure Storage integration

In order to perform this test, you need an Azure Storage account
with read/write access.
Both the storage account and associated key or SAS token (without
leading questionmark) need to be
passed to snakemake at runtime, by exporting
environment variables `AZURE_ACCOUNT` and either `AZURE_KEY` or
`SAS_TOKEN`.

Furthermore, in the storage account, a container "snakemake-test"
needs to be created prior to running the test.

And lastly, a local file called `test.txt.gz` needs to be created.

