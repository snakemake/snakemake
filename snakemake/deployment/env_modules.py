__author__ = "Johannes Köster"
__copyright__ = "Copyright 2015-2019, Johannes Köster"
__email__ = "johannes.koester@uni-due.de"
__license__ = "MIT"


class EnvModules:
    def __init__(self, *module_names):
        self.names = module_names

    def shellcmd(self, cmd):
        """Return shell command with given modules to load, whilst
           ensuring a clean environment.
        """
        to_load = " ".join(name for name in self.names)
        return "module purge && module load {}".format(to_load)

    def __str__(self):
        return ", ".join(self.names)
