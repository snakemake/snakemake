__author__ = "Sebastian Kurscheid (sebastian.kurscheid@anu.edu.au)"
__license__ = "MIT"
__date__ = "2018-09-015"

from snakemake.exceptions import MissingInputException
import os

from snakemake.remote.AzureStorage import RemoteProvider as AzureRemoteProvider
account_key='khTRU9zcTi6lTrBp7IsswmZRGrzEIOxBXiF1AI/uCAdbFK0kgZYmz6W/IETswtSZg1EFNOoNhYL75R2Ir1bEMQ=='
account_name='cellcycledata'

AS = AzureRemoteProvider(account_name=account_name, account_key=account_key)
file_names = AS.glob_wildcards("experiment/{fn}.gz")

rule glob_test:
    input:
        expand("local_data/{fn}.gz", fn=file_names))

rule make_data:
    input:
        AS.remote("experiment/aks_test/{fn}.gz")
    output:
        "local_data/{fn}.gz"
    shell:
        "mv {input[0]} {output[0]}"

rule all:
    input:
        AS.remote("e-mtab-6967/E-MTAB-6967.sdrf.txt")
    run:
        outputName = os.path.basename(input[0])
        shell("cp {input} {outputName}")