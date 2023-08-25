from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

from snakemake.logging import logger


@dataclass
class WorkdirHandler:
    workdir: Optional[Path] = None
    olddir: Optional[Path] = field(init=False, default=None)

    def change_to(self):
        if self.workdir is not None:
            self.olddir = Path.cwd()
            if not self.workdir.exists():
                logger.info(f"Creating specified working directory {self.workdir}.")
                self.workdir.mkdir(parents=True)
            self.workdir.chdir()

    def change_back(self):
        if self.workdir is not None:
            self.olddir.chdir()