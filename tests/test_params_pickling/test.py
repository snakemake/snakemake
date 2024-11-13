import numpy as np
import polars as pl
import pandas as pd

assert isinstance(snakemake.params.testnp, np.ndarray)
assert isinstance(snakemake.params.testpl, pl.DataFrame)
assert isinstance(snakemake.params.testpd, pd.DataFrame)

snakemake.params.testpd.to_csv(snakemake.output.pd, sep="\t")
snakemake.params.testpl.write_csv(snakemake.output.pl, separator="\t")
np.savetxt(snakemake.output.np, snakemake.params.testnp, delimiter="\t")