__author__ = "Johannes Köster"
__copyright__ = "Copyright 2015-2019, Johannes Köster"
__email__ = "johannes.koester@uni-due.de"
__license__ = "MIT"


class EnvModules:
    def __init__(self, *module_names):
        self.names = module_names

    def shellcmd(self, cmd):
        """Return shell command with given modules loaded."""
        loads = " ".join(self._load_module(name) for name in self.names)
        return "{} {}".format(loads, cmd)

    def _load_module(self, name):
        return "module load {};".format(name)

    def __str__(self):
        return ", ".join(self.names)
