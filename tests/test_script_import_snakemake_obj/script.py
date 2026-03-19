from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from snakemake.io.container import snakemake


with open(snakemake.output[0], "w") as f:
    f.write("Hello world!")