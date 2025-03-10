from snakemake.shell import shell
print("test", file=open(snakemake.output[0], "w"))
shell("echo test")