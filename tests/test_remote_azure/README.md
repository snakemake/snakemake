# Instruction for testing of Azure Storage integration
* in order to perform this test, an Azure Storage Account is required
* Both the storage account and associated key need to be passed to snakemake at runtime
* currently this is solved by setting and exporting environment variables called
** $AZURE_ACCOUNT
** $AZURE_KEY
* furthermore, in the storage account, a container "snakemake-test" needs to be created prio to running the test

