from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

from snakemake.logging import logger


@dataclass
class WorkdirHandler:
    workdir: Optional[Path] = None
    olddir: Path = field(init=False)

    def __enter__(self):
        if self.workdir is not None:
            self.olddir = Path.cwd()
            if not self.workdir.exists():
                logger.info(f"Creating specified working directory {self.workdir}.")
                self.workdir.mkdir(parents=True)
            self.workdir.chdir()

    def __exit__(self, exc_type, exc_value, traceback):
        if self.workdir is not None:
            self.olddir.chdir()