# Snakemake Tutorial with ORAS

This is a derivative of the [Snakemake Tutorial](https://snakemake.readthedocs.io/en/stable/tutorial/tutorial.html)
that combines the [data needed](https://github.com/snakemake/snakemake-tutorial-data/tree/master) with an OCI registry,
and downloads using ORAS (OCI Registry as Storage). Given that you have a conda environment ready, you can run as follows:

```bash
# if you have mamba
$ snakemake --cores 1 --use-conda

# if you have conda
$ snakemake --cores 1 --use-conda --conda-frontend=conda
```

If you don't have conda or mamba (and remove the flag) you'll need samtools, bwa, and the other dependency binaries
available locally. You can also run this example in a container:

```bash
$ docker run -v $PWD:/code -it snakemake/snakemake 
```
And then (from the same example directory):

```bash
$ cd /code
$ pip install oras
$ snakemake --cores 1 --use-conda
```