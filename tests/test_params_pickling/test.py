import numpy as np
import polars as pl
import pandas as pd

assert type(snakemake.params.testnp) == np.ndarray
assert type(snakemake.params.testpl) == pl.DataFrame
assert type(snakemake.params.testpd) == pd.DataFrame

snakemake.params.testpd.to_csv(snakemake.output.pd, sep="\t")
snakemake.params.testpl.write_csv(snakemake.output.pl, separator="\t")
np.savetxt(snakemake.output.np, snakemake.params.testnp, delimiter="\t")