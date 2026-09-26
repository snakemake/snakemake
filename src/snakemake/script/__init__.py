"""Runtime context for Snakemake scripts."""

from snakemake.iocontainers import Snakemake

is_script = False
snakemake: Snakemake
