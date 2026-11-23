import platform
import PIL
from snakemake.utils import makedirs

assert snakemake.params.x == 1

with open("version.txt", "w") as f:
    f.write(platform.python_version())
