__author__ = "Johannes Köster"
__copyright__ = "Copyright 2022, Johannes Köster"
__email__ = "johannes.koester@uni-due.de"
__license__ = "MIT"


import hashlib


class EnvModules:
    def __init__(self, *module_names):
        self.names = module_names

    def shellcmd(self, cmd):
        """Return shell command with given modules loaded."""
        return "module purge && module load {to_load}; {cmd}".format(
            to_load=" ".join(self.names), cmd=cmd
        )

    def __str__(self):
        return ", ".join(self.names)

    @property
    def hash(self) -> str:
        md5hash = hashlib.md5()
        for name in self.names:
            md5hash.update(name.encode())
        return md5hash.hexdigest()
